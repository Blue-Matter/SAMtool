# Rapid Conditioning Model (RCM)

Intended for conditioning operating models for MSEtool. For data-limited
stocks, this function can generate a range of potential depletion
scenarios inferred from sparse data. From a historical time series of
total catch or effort, and potentially age/length compositions and
multiple indices of abundance, the RCM returns a range of values for
depletion, selectivity, unfished recruitment (R0), historical fishing
effort, and recruitment deviations for the operating model. This is done
by sampling life history parameters provided by the user and fitting a
statistical catch-at-age model (with the predicted catch equal to the
observed catch). Alternatively one can do a single model fit and sample
the covariance matrix to generate an operating model with uncertainty
based on the model fit. Either a full catch (conditioned on catch) or
effort (conditioned on effort) time series is needed but missing data
(as NAs) are allowed for all other data types. `check_RCMdata` evaluates
whether the inputs in the S4 RCMdata object are correctly formatted.

## Usage

``` r
check_RCMdata(RCMdata, OM, condition = "catch", silent = FALSE)

RCM(OM, data, ...)

# S4 method for class 'OM,RCMdata'
RCM(
  OM,
  data,
  condition = "catch",
  selectivity = "logistic",
  s_selectivity = NULL,
  LWT = list(),
  comp_like = c("multinomial", "lognormal", "mvlogistic", "dirmult1", "dirmult2"),
  prior = list(),
  max_F = 3,
  cores = 1L,
  integrate = FALSE,
  mean_fit = FALSE,
  drop_nonconv = FALSE,
  drop_highF = FALSE,
  control = list(iter.max = 2e+05, eval.max = 4e+05),
  start = list(),
  map = list(),
  silent = FALSE,
  ...
)

# S4 method for class 'OM,list'
RCM(
  OM,
  data,
  condition = "catch",
  selectivity = "logistic",
  s_selectivity = NULL,
  LWT = list(),
  comp_like = c("multinomial", "lognormal", "mvlogistic", "dirmult1", "dirmult2"),
  ESS = c(30, 30),
  prior = list(),
  max_F = 3,
  cores = 1L,
  integrate = FALSE,
  mean_fit = FALSE,
  drop_nonconv = FALSE,
  drop_highF = FALSE,
  control = list(iter.max = 2e+05, eval.max = 4e+05),
  start = list(),
  map = list(),
  silent = FALSE,
  ...
)

# S4 method for class 'OM,Data'
RCM(
  OM,
  data,
  condition = "catch",
  selectivity = "logistic",
  s_selectivity = NULL,
  LWT = list(),
  comp_like = c("multinomial", "lognormal", "mvlogistic", "dirmult1", "dirmult2"),
  ESS = c(30, 30),
  prior = list(),
  max_F = 3,
  cores = 1L,
  integrate = FALSE,
  mean_fit = FALSE,
  drop_nonconv = FALSE,
  drop_highF = FALSE,
  control = list(iter.max = 2e+05, eval.max = 4e+05),
  start = list(),
  map = list(),
  silent = FALSE,
  ...
)

# S4 method for class 'list,RCMdata'
RCM(
  OM,
  data,
  condition = "catch",
  selectivity = "logistic",
  s_selectivity = NULL,
  LWT = list(),
  comp_like = c("multinomial", "lognormal", "mvlogistic", "dirmult1", "dirmult2"),
  prior = list(),
  max_F = 3,
  integrate = FALSE,
  control = list(iter.max = 2e+05, eval.max = 4e+05),
  start = list(),
  map = list(),
  silent = FALSE,
  ...
)
```

## Arguments

- RCMdata:

  An [RCMdata](https://samtool.openmse.com/reference/RCMdata-class.md)
  object.

- OM:

  An object of class
  [MSEtool::OM](https://msetool.openmse.com/reference/OM-class.html)
  that specifies natural mortality (M), growth (Linf, K, t0, a, b),
  stock-recruitment relationship, steepness, maturity parameters (L50
  and L50_95), and standard deviation of recruitment variability (Perr).
  Alternatively, provide a named list of biological inputs, see
  "StockPars" section below.

- condition:

  String to indicate whether the RCM is conditioned on "catch" (where F
  are estimated parameters), "catch2" (where F is solved internally
  using Newton's method), or "effort" (F is proportional to an index
  series in `data@Ehist`. Can be fleet specific, in which case some
  combination of "catch" and "effort" are permissible.

- silent:

  Logical to indicate whether informative messages will be reported to
  console.

- data:

  Data inputs formatted in a
  [RCMdata](https://samtool.openmse.com/reference/RCMdata-class.md)
  (preferred) or
  [MSEtool::Data](https://msetool.openmse.com/reference/Data-class.html)
  object. Use of a list is deprecated. See Data section below.

- ...:

  Other arguments to pass in for starting values of parameters and
  fixing parameters. See details.

- selectivity:

  A character vector of length nfleet to indicate `"logistic_length"`,
  `"dome_length"`, `"logistic_age"`, `"dome_age"`, or `"free"`
  selectivity for each fleet in `Chist`. If there is time-varying
  selectivity, this is a character vector of length nsel_block (see Data
  section below). "free" indicates independent selectivity parameters
  for each age, and additional modifications for fixing selectivity
  parameters will likely be needed. See Additional arguments section.

- s_selectivity:

  A vector of length nsurvey to indicate the selectivity of the
  corresponding columns in `data$Index`. Use `"B"` for total biomass, or
  `"SSB"` for spawning biomass (by default, "B" is used). Use numbers if
  the survey selectivity follows a fleet (corresponding to the columns
  in data\$Chist, e.g., 1 = first fleet/column and so on). If the survey
  selectivity is otherwise independent of anything else in the model,
  use `"logistic_length"`, `"dome_length"`, `"logistic_age"`,
  `"dome_age"`, or `"free"` to specify the functional form of
  selectivity, and see Additional arguments section for setup of survey
  selectivity parameters and Articles section for more information.

- LWT:

  A named list of likelihood weights for the RCM. See below.

- comp_like:

  A string indicating the statistical distribution for the composition
  data, either `"multinomial"` (default), `"lognormal"`, `"mvlogistic"`
  (multivariate logistic), `"dirmult1"` (Dirichlet multinomial, linear
  version), or `"dirmult2"` (saturating version; see Thorson et al.
  2017).

- prior:

  A named list for the parameters of any priors to be added to the
  model. See below.

- max_F:

  The maximum F for any fleet in the scoping model (higher F's in the
  model are penalized in the objective function). This argument will
  also update `OM@maxF`. See also `drop_highF`.

- cores:

  Integer for the number of CPU cores (set greater than 1 for parallel
  processing).

- integrate:

  Logical, whether to treat recruitment deviations as penalized
  parameters in the likelihood (FALSE) or random effects to be
  marginalized out of the likelihood (TRUE).

- mean_fit:

  Logical, whether to run an additional with mean values of life history
  parameters from the OM.

- drop_nonconv:

  Logical, whether to drop non-converged fits of the RCM, including fits
  where F = NA.

- drop_highF:

  Logical, whether to drop fits of the RCM where F = `max_F`.

- control:

  A named list of arguments (e.g, max. iterations, etc.) for
  optimization, to be passed to the control argument of
  [`stats::nlminb()`](https://rdrr.io/r/stats/nlminb.html).

- start:

  A list of starting values for the TMB model. See details.

- map:

  A list of `map` argument to TMB models to override defaults. See
  [MakeADFun](https://rdrr.io/pkg/TMB/man/MakeADFun.html) and details.

- ESS:

  A vector of length two. A shortcut method to setting the maximum
  multinomial sample size of the age and length compositions. Not used
  when data are provided in a
  [RCMdata](https://samtool.openmse.com/reference/RCMdata-class.md)
  object.

## Value

An object of class
[RCModel](https://samtool.openmse.com/reference/RCModel-class.md) (see
link for description of output).

`check_RCMdata` returns a list of updated RCMdata object, OM, and
StockPars and FleetPars from the Hist object generated from the OM.

## Details

Fleet selectivity is fixed to values sampled from `OM` if no age or
length compositions are provided.

Survey selectivity is estimable only if `IAA` or `IAL` is provided.
Otherwise, the selectivity should be mirrored to a fleet (vulnerable
biomass selectivity) or indexed to total or spawning biomass (see
`s_selectivity`).

Parameters that were used in the fitting model are placed in the
`RCM@OM@cpars` list.

If the operating model `OM` uses time-varying growth or M, then those
trends will be used in the RCM as well. Non-stationary productivity
creates ambiguity in the calculation and interpretation of depletion and
MSY reference points.

The easiest way to turn off time-varying growth/M is by setting:
`OM@Msd <- OM@Linfsd <- OM@Ksd <- c(0, 0)`.

To play with alternative fits by excluding indices, for example, or
other optional data, set the corresponding likelihood weight to zero.
The model will still generate the inferred index but the data won't
enter the likelihood. See section on likelihood weights.

## Online Documentation

Several articles are available for RCM:

- [General overview of approach](https://openmse.com/tutorial-rcm/)

- [Mathematical description](https://openmse.com/tutorial-rcm-eq/)

- [Setup of selectivity settings and index
  catchability](https://openmse.com/tutorial-rcm-select/) (useful for
  more data-rich cases)

- [Description of
  priors](https://openmse.com/features-assessment-models/5-priors/)

## Priors

The following priors can be added as a named list, e.g.,
`prior = list(M = c(0.25, 0.15), h = c(0.7, 0.1)`. For each parameter
below, provide a vector of values as described:

- `R0`:

  A vector of length 3. The first value indicates the distribution of
  the prior: `1` for lognormal, `2` for uniform on `log(R0)`, `3` for
  uniform on R0. If lognormal, the second and third values are the prior
  mean (in normal space) and SD (in log space). Otherwise, the second
  and third values are the lower and upper bounds of the uniform
  distribution (values in normal space).

- `h`:

  A vector of length 2 for the prior mean and SD, both in normal space.
  Beverton-Holt steepness uses a beta distribution, while Ricker
  steepness uses a normal distribution.

- `M`:

  A vector of length 2 for the prior mean (in normal space) and SD (in
  log space). Lognormal prior.

- `q`:

  A matrix for nsurvey rows and 2 columns. The first column is the prior
  mean (in normal space) and the second column for the SD (in log
  space). Use `NA` in rows corresponding to indices without priors.

See online documentation for more details.

## Data

One of indices, age compositions, or length compositions should be
provided in addition to the historical catch or effort. Not all
arguments are needed to run the model (some have defaults, while others
are ignored if not applicable depending on the data provided).

The `data` variable can be an object of class
[RCMdata](https://samtool.openmse.com/reference/RCMdata-class.md). See
help file for description of inputs.

Alternatively, the `data` input can be a
[MSEtool::Data](https://msetool.openmse.com/reference/Data-class.html)
S4 object which will retrieve data from the following slots:

- `Data@Cat`:

  catch series (single fleet with the Data S4 object)

- `Data@Effort`:

  effort series

- `Data@CAA`:

  fishery age composition

- `Data@CAL`, `Data@CAL_mids`:

  fishery length composition and corresponding length bins

- `Data@Ind`, `Data@SpInd`, `Data@VInd`, `Data@AddInd`:

  indices of abundance

- `Data@CV_Ind`, `Data@CV_SpInd`, `Data@CV_VInd`, `Data@CV_AddInd`:

  annual coefficients of variation for the corresponding indices of
  abundance. CVs will be converted to lognormal standard deviations.

- `Data@ML`:

  fishery mean lengths

- `Data@AddIndV`, `Data@AddIndType`, `Data@AddIunits`:

  Additional information for indices in `Data@AddInd`: selectivity and
  units (i.e., biomass or abundance).

There is no slot in the Data S4 object for the equilibrium catch/effort.
These can be passed directly in the function call, i.e.,
`RCM(OM, Data, C_eq = C_eq, ...)`.

## Data list (deprecated)

Use of a list is deprecated. For backwards compatibility, here is the
list of supported entries:

- `Chist`:

  A vector of historical catch, should be of length OM@nyears. If there
  are multiple fleets: a matrix of `OM@nyears` rows and `nfleet`
  columns. Ideally, the first year of the catch series represents
  unfished conditions (see also `C_eq`).

- `C_sd`:

  A vector or matrix of standard deviations (lognormal distribution) for
  the catches in `Chist`. If not provided, the default is 0.01. Only
  used if `condition = "catch"`.

- `Ehist`:

  A vector of historical effort, should be of length `OM@nyears` (see
  also `E_eq`).

- `Index`:

  A vector of values of an index (of length `OM@nyears`). If there are
  multiple indices: a matrix of historical indices of abundances, with
  rows indexing years and columns indexing the index.

- `I_sd`:

  A vector or matrix of standard deviations (lognormal distribution) for
  the indices corresponding to the entries in `Index`. If not provided,
  this function will use values from `OM@Iobs`.

- `I_type`:

  Obsolete as of version 2.0. See `s_selectivity` argument.

- `CAA`:

  Fishery age composition matrix with `nyears` rows and `OM@maxage+1`
  columns. If multiple fleets: an array with dimension:
  `nyears, OM@maxage, and nfleet`.

- `CAL`:

  Fishery length composition matrix with nyears rows and columns
  indexing the length bin. If multiple fleets: an array with dimension:
  `nyears, length bins, and nfleet`.

- `MS`:

  A vector of fishery mean size (MS, either mean length or mean weight)
  observations (length `OM@nyears`), or if multiple fleets: matrix of
  dimension: `nyears, nfleet`. Generally, mean lengths should not be
  used if `CAL` is also provided, unless mean length and length comps
  are independently sampled.

- `MS_type`:

  A character (either `"length"` (default) or `"weight"`) to denote the
  type of mean size data.

- `MS_cv`:

  The coefficient of variation of the observed mean size. If there are
  multiple fleets, a vector of length `nfleet`. Default is 0.2.

- `s_CAA`:

  Survey age composition data, an array of dimension
  `nyears, maxage+1, nsurvey`.

- `s_CAL`:

  Survey length composition data, an array of dimension
  `nyears, length(length_bin), nsurvey`.

- `length_bin`:

  A vector for the midpoints of the length bins for `CAL` and `s_CAL`.
  All bin widths should be equal in size.

- `C_eq`:

  A numeric vector of length `nfleet` for the equilibrium catch for each
  fleet in `Chist` prior to the first year of the operating model. Zero
  (default) implies unfished conditions in year one. Otherwise, this is
  used to estimate depletion in the first year of the data.
  Alternatively, if one has a full CAA matrix, one could instead
  estimate "artificial" rec devs to generate the initial numbers-at-age
  (and hence initial depletion) in the first year of the model (see
  additional arguments).

- `C_eq_sd`:

  A vector of standard deviations (lognormal distribution) for the
  equilibrium catches in `C_eq`. If not provided, the default is 0.01.
  Only used if `condition = "catch"`.

- `E_eq`:

  The equilibrium effort for each fleet in `Ehist` prior to the first
  year of the operating model. Zero (default) implies unfished
  conditions in year one. Otherwise, this is used to estimate depletion
  in the first year of the data.

- `abs_I`:

  Optional, an integer vector to indicate which indices are in absolute
  magnitude. Use 1 to set `q = 1`, otherwise use 0 to estimate q.

- `I_units`:

  Optional, an integer vector to indicate whether indices are biomass
  based (1) or abundance-based (0). By default, all are biomass-based.

- `age_error`:

  Optional, a square matrix of maxage + 1 rows and columns to specify
  ageing error. The aa-th column assigns a proportion of the true age in
  the a-th row to observed age. Thus, all rows should sum to 1. Default
  is an identity matrix (no ageing error).

- `sel_block`:

  Optional, for time-varying fleet selectivity (in time blocks), a
  integer matrix of `nyears` rows and `nfleet` columns to assigns a
  selectivity function to a fleet for certain years.

## StockPars

When an operating model is provided, the RCM function will generally fit
to each simulation of biological parameters.

Alternatively for a single fit to data independent of any operating
model, provide a named list containing the following (naming conventions
follow internal operating model variables):

- `SRrel` Integer, stock-recruit function (1 = Beverton-Holt, 2 =
  Ricker, 3 = Mesnil-Rochet hockey stick)

- `R0` Numeric, starting value for unfished recruitment parameter

- `M_ageArray` Matrix `[maxage+1, nyears]` for natural mortality

- `Len_age` Matrix `[maxage+1, nyears + 1]` for length at age

- `Linf` Numeric. Asymptotic length. Only used for the upper bound for
  the size of full selectivity (if selectivity functions are
  length-based)

- `LatASD` Matrix `[maxage+1, nyears + 1]` for the standard deviation in
  length at age

- `Wt_age` Matrix `[maxage+1, nyears + 1]` for stock weight at age

- `Mat_age` Matrix `[maxage+1, nyears + 1]` for maturity at age

- `Fec_Age` Matrix `[maxage+1, nyears + 1]` for fecundity at age.
  Frequently the product of maturity and weight at age

- `ageMarray` Numeric, age of 50 percent maturity. Used to average the
  initial years for the unfished replacement line of the stock recruit
  relationship and steepness/R0. Irrelevant if fecundity and natural
  mortality are not time-varying (set to 1).

- `spawn_time_frac` Numeric, fraction of the year when spawning occurs

- `hs` Numeric, steepness of the stock recruit relationship

- `procsd` Numeric, lognormal recruitment deviation standard deviation

## Additional arguments

For `RCM`, additional arguments can be passed to the model via `...`:

- `plusgroup`:

  Logical for whether the maximum age is a plusgroup or not. By default,
  TRUE.

- `fix_dome`:

  Logical for whether the dome selectivity parameter for fleets is
  fixed. Used primarily for backwards compatibility, this is overridden
  by the `map` argument.

- `resample`:

  Logical, whether the OM conditioning parameters (recruitment, fishing
  mortality, SSB, selectivity, etc.) are obtained by sampling the
  Hessian matrix from a single model fit. By default FALSE. This feature
  requires identical biological parameters among simulations.

- `pbc_recdev`:

  Vector of length nyears. Proportion of the bias correction to apply
  annually to the recruitment deviations (if estimated). The bias
  correction from logspace to normal space is
  `exp(log_rec_dev[y] - 0.5 * pbc_recdev[y] * sigmaR^2)`. Default
  proportion is 1.

- `pbc_earlyrecdev`:

  Vector of length maxage. Proportion of the bias correction to apply to
  the abundance deviations in the first year of the model (if
  estimated). The bias correction from logspace to normal space is
  `exp(log_early_rec_dev[a] - 0.5 * pbc_recdev[a] * sigmaR^2)`. Default
  proportion is 1.

## start list

Starting values can be specified in a named list for the following:

- `vul_par`:

  A matrix of 3 rows and nfleet columns for starting values for fleet
  selectivity. The three rows correspond to LFS (length of full
  selectivity), L5 (length of 5 percent selectivity), and Vmaxlen
  (selectivity at length Linf). By default, the starting values are
  values from the OM object. If any `selectivity = "free"`, then this
  matrix needs to be of `maxage+1` rows where the row specifies the
  selectivity at age. See Articles section.

- `ivul_par`:

  A matrix of 3 rows and nsurvey columns for starting values for fleet
  selectivity. Same setup as `vul_par`. Values in the column are ignored
  if `s_selectivity` is mapped to a fishing fleet (add NA placeholders
  in that case). If any `s_selectivity = "free"`, then this matrix needs
  to be of `maxage+1` rows where the row specifies the selectivity at
  age.

- `log_rec_dev`:

  A numeric vector of length `nyears` for the starting values of the
  log-recruitment deviations.

- `log_early_rec_dev`:

  A numeric vector of length `OM@maxage` for the starting values of the
  recruitment deviations controlling the abundance-at-age in the first
  year of the model.

- `q`:

  A numeric vector of length nsurvey for index catchability. See [online
  article](https://openmse.com/tutorial-rcm-select/) for more
  information.

## map list

Parameters can be fixed with the map argument (also a named list,
corresponding to the start list). Each vector or matrix in the map
argument will be the same dimension as in the start entry. If an entry
is `NA`, the corresponding parameter is fixed in the model to the
starting value. Otherwise, an integer for each independent parameter,
i.e., shared or mirrored parameters get the same integer entry.

- `vul_par`:

  An integer matrix of the same dimension as `start$vul_par`. By
  default, selectivity is fixed if there are no age or length
  composition for that fleet or survey, otherwise estimated. Unused
  cells in the `start$vul_par` matrix should be given NA in the map
  matrix.

- `ivul_par`:

  The map argument for the survey selectivity parameters (same dimension
  as `start$ivul_par`). Placeholder parameters should have a map value
  of NA.

- `log_early_rec_dev`:

  A vector of length `OM@maxage` that indexes which recruitment deviates
  for the cohorts in the first year of the model are fixed (using NA) or
  estimated (a separate integer). By default, no deviates are estimated
  (all are NA).

- `log_rec_dev`:

  A vector of length `OM@nyears` that indexes which recruitment deviates
  are fixed (using NA) or estimated (a separate integer). By default,
  all these deviates are estimated.

- `q`:

  A vector of length `nsurvey` for index catchability. q should be an
  estimated parameter when sharing across surveys (perhaps with
  differing selectivity). Otherwise, it is solved analytically where
  individual parameters are independent of other indices. Use
  `RCMdata@abs_I` for fixing the catchability to 1. See [online
  article](https://openmse.com/tutorial-rcm-select/) for more
  information.

## Likelihood weights

`LWT` is an optional named list containing the likelihood weights
(values \>= 0) with the possible options:

- `Chist, CAA, CAL, MS, C_eq`: A vector of length nfleet for each.

- `Index, IAA, IAL`: A vector of length nsurvey for each.

By default, all likelihood weights are equal to one if not specified by
the user.

Annual multinomial sample sizes for the age and length comps can now be
provided directly in the
[RCMdata](https://samtool.openmse.com/reference/RCMdata-class.md)
object. For a list or
[MSEtool::Data](https://msetool.openmse.com/reference/Data-class.html)
object, use the `ESS` argument.

## References

Thorson et al. 2017. Model-based estimates of effective sample size in
stock assessment models using the Dirichlet-multinomial distribution.
Fish. Res. 192:84-93.
[doi:10.1016/j.fishres.2016.06.005](https://doi.org/10.1016/j.fishres.2016.06.005)

## See also

[plot.RCModel](https://samtool.openmse.com/reference/plot.RCModel.md)
[RCModel](https://samtool.openmse.com/reference/RCModel-class.md)
[compare_RCM](https://samtool.openmse.com/reference/plot.RCModel.md)
[pcod](https://samtool.openmse.com/reference/pcod.md)
[RCM2MOM](https://samtool.openmse.com/reference/RCM2MOM.md)
[posterior](https://samtool.openmse.com/reference/posterior.md)

## Author

Q. Huynh

## Examples

``` r
# \donttest{ 
# An example that conditions a Pacific cod operating model. There are 48 simulations, 
# where values of natural mortality and steepness are sampled from distributions. 
# The model is fitted with priors on the index catchability. Maturity and selectivity 
# are knife-edge at the age of 2 years. See online tutorial for more information.

data(pcod) 
mat_ogive <- pcod$OM@cpars$Mat_age[1, , 1]
out <- RCM(OM = pcod$OM, data = pcod$data, 
           condition = "catch", mean_fit = TRUE,
           selectivity = "free", s_selectivity = rep("SSB", ncol(pcod$data@Index)),
           start = list(vul_par = matrix(mat_ogive, length(mat_ogive), 1)),
           map = list(vul_par = matrix(NA, length(mat_ogive), 1),
                      log_early_rec_dev = rep(1, pcod$OM@maxage)),
           prior = pcod$prior)
#> ℹ Checking data...
#> ✔ 1 fleet(s) detected.
#> ✔ RCM is conditioned on:
#> ✔ Fleet 1: catch
#> ✔ 65 years of data detected.
#> ✔ First year in model: 1956
#> ✔ Last year in model: 2020
#> ✔ 5 survey(s) detected.
#> ✔ Checking OM and getting biological parameters...
#> ✔ Mean weight data found.
#> ✔ Passing user arguments (LWT, map, start, prior, etc.) to RCMdata@Misc..
#> ✔ Maximum F in RCM will be 3. OM@maxF is also updated.
#>   
#> ℹ No fishery length or age compositions were provided. Selectivity is fixed to values from OM.
#>   
#> ℹ Fishery selectivity setup:
#> ℹ Fleet 1: individual parameters at age (free)
#>   
#> ℹ Index selectivity setup:
#> ℹ Index 1: spawning biomass
#> ℹ Index 2: spawning biomass
#> ℹ Index 3: spawning biomass
#> ℹ Index 4: spawning biomass
#> ℹ Index 5: spawning biomass
#>   
#> ✔ Beverton-Holt stock-recruitment relationship used.
#> ✔ Prior for q found.
#> ℹ Fitting model (48 simulations) ...
#> ℹ Generating additional model fit from mean values of parameters in the operating model...
#> ✔ Updating operating model:
#>   
#> ✔ Range of unfished age-0 recruitment (OM@cpars$R0): 6383.47 - 13782.79
#> ✔ Range of initial spawning depletion: 0.44 - 1.45
#> ✔ Range of spawning depletion (OM@cpars$D): 0.14 - 0.45
#> ✔ Historical F set with OM@cpars$Find and OM@cpars$qs.
#> ✔ Annual selectivity at age set in OM@cpars$V. Projection period uses selectivity of last historical year.
#> ✔ RCMdata length bins will be added to OM.
#> ✔ Recruitment standard deviation set in OM@cpars$Perr: 0.8 - 0.8
#> ✔ Historical recruitment deviations set in OM@cpars$Perr_y.
#> ✔ Range of recruitment autocorrelation OM@AC: 0.21 - 0.31
#> ✔ Future recruitment deviations in OM@cpars$Perr_y sampled with autocorrelation.
#> ✔ Growth, maturity, natural mortality, and stock recruit parameters from RCM are set in OM@cpars.
#>   
#> ℹ Adding some RCMdata inputs into OM@cpars$Data:
#>   
#> ✔ Historical catch data added to OM@cpars$Data@Cat.
#> ✔ Historical indices added to OM@cpars$Data@AddInd.
#> ✔ Complete.
plot(out, s_name = colnames(pcod$data@Index))
#> ℹ Rendering markdown file: /tmp/RtmpNbrYWx/RCM.Rmd
#> ℹ See help(plot.RCModel) to adjust report and file directory.
#> ✔ Rendered file: /tmp/RtmpNbrYWx/RCM.html

# Alternative OM with age-3 maturity and selectivity instead.
out_age3 <- local({
  pcod$OM@cpars$Mat_age[, 2, ] <- 0
  mat_ogive_age3 <- pcod$OM@cpars$Mat_age[1, , 1]
  RCM(OM = pcod$OM, data = pcod$data, 
      condition = "catch", mean_fit = TRUE,
      selectivity = "free", s_selectivity = rep("SSB", ncol(pcod$data@Index)),
      start = list(vul_par = matrix(mat_ogive_age3, length(mat_ogive_age3), 1)),
      map = list(vul_par = matrix(NA, length(mat_ogive_age3), 1),   
                 log_early_rec_dev = rep(1, pcod$OM@maxage)),
      prior = pcod$prior)
})
#> ℹ Checking data...
#> ✔ 1 fleet(s) detected.
#> ✔ RCM is conditioned on:
#> ✔ Fleet 1: catch
#> ✔ 65 years of data detected.
#> ✔ First year in model: 1956
#> ✔ Last year in model: 2020
#> ✔ 5 survey(s) detected.
#> ✔ Checking OM and getting biological parameters...
#> ✔ Mean weight data found.
#> ✔ Passing user arguments (LWT, map, start, prior, etc.) to RCMdata@Misc..
#> ✔ Maximum F in RCM will be 3. OM@maxF is also updated.
#>   
#> ℹ No fishery length or age compositions were provided. Selectivity is fixed to values from OM.
#>   
#> ℹ Fishery selectivity setup:
#> ℹ Fleet 1: individual parameters at age (free)
#>   
#> ℹ Index selectivity setup:
#> ℹ Index 1: spawning biomass
#> ℹ Index 2: spawning biomass
#> ℹ Index 3: spawning biomass
#> ℹ Index 4: spawning biomass
#> ℹ Index 5: spawning biomass
#>   
#> ✔ Beverton-Holt stock-recruitment relationship used.
#> ✔ Prior for q found.
#> ℹ Fitting model (48 simulations) ...
#> ℹ Generating additional model fit from mean values of parameters in the operating model...
#> ✔ Updating operating model:
#>   
#> ✔ Range of unfished age-0 recruitment (OM@cpars$R0): 6383.47 - 13782.79
#> ✔ Range of initial spawning depletion: 0.44 - 1.45
#> ✔ Range of spawning depletion (OM@cpars$D): 0.14 - 0.45
#> ✔ Historical F set with OM@cpars$Find and OM@cpars$qs.
#> ✔ Annual selectivity at age set in OM@cpars$V. Projection period uses selectivity of last historical year.
#> ✔ RCMdata length bins will be added to OM.
#> ✔ Recruitment standard deviation set in OM@cpars$Perr: 0.8 - 0.8
#> ✔ Historical recruitment deviations set in OM@cpars$Perr_y.
#> ✔ Range of recruitment autocorrelation OM@AC: 0.21 - 0.31
#> ✔ Future recruitment deviations in OM@cpars$Perr_y sampled with autocorrelation.
#> ✔ Growth, maturity, natural mortality, and stock recruit parameters from RCM are set in OM@cpars.
#>   
#> ℹ Adding some RCMdata inputs into OM@cpars$Data:
#>   
#> ✔ Historical catch data added to OM@cpars$Data@Cat.
#> ✔ Historical indices added to OM@cpars$Data@AddInd.
#> ✔ Complete.
  
compare_RCM(out, out_age3, scenario = list(names = c("Age-2 maturity", "Age-3 maturity")),
            s_name = colnames(pcod$data@Index))
#> ℹ Rendering markdown file: /tmp/RtmpNbrYWx/compare_RCM.Rmd
#> ✔ Rendered file: /tmp/RtmpNbrYWx/compare_RCM.html
             
Hist <- runMSE(out@OM, Hist = TRUE)            
#> → Loading operating model
#> ✔ Note: Maximum age (10) is lower than assuming 1% of cohort survives to maximum age (24)
#> → Optimizing for user-specified movement
#> → Calculating MSY reference points for each year
#> → Calculating historical stock and fishing dynamics
#> → Calculating per-recruit reference points
#> → Calculating B-low reference points
#> → Calculating reference yield - best fixed F strategy
#> → Simulating observed data
#> → Updating Simulated Data with Real Data from `OM@cpars$Data`
#> Warning: An error occurred in calculating statistical properties of fit to Additional Index 3 (possibly because there was only one observed data point). 
#> Using the index observation error for slot `Ind` from `Obs` object (or possibly conditioned if `cpars$Data@Ind` was provided).
#> Use `cpars$AddIerr` to manually set the observation error.
#> → Returning historical simulations
# } 
```
