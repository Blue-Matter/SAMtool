# Class-`Assessment`

An S4 class that contains assessment output. Created from a function of
class `Assess`.

## Slots

- `Model`:

  Name of the assessment model.

- `Name`:

  Name of Data object.

- `conv`:

  Logical. Whether the assessment model converged (defined by whether
  TMB returned a positive-definite covariance matrix for the model).

- `UMSY`:

  Estimate of exploitation at maximum sustainable yield.

- `FMSY`:

  Estimate of instantaneous fishing mortality rate at maximum
  sustainable yield.

- `MSY`:

  Estimate of maximum sustainable yield.

- `BMSY`:

  Biomass at maximum sustainable yield.

- `SSBMSY`:

  Spawning stock biomass at maximum sustainable yield.

- `VBMSY`:

  Vulnerable biomass at maximum sustainable yield.

- `B0`:

  Biomass at unfished equilibrium.

- `R0`:

  Recruitment at unfished equilibrium.

- `N0`:

  Abundance at unfished equilibrium.

- `SSB0`:

  Spawning stock biomass at unfished equilibrium.

- `VB0`:

  Vulnerable biomass at unfished equilibrium.

- `h`:

  Steepness.

- `U`:

  Time series of exploitation.

- `U_UMSY`:

  Time series of relative exploitation.

- `FMort`:

  Time series of instantaneous fishing mortality.

- `F_FMSY`:

  Time series of fishing mortality relative to MSY.

- `B`:

  Time series of biomass.

- `B_BMSY`:

  Time series of biomass relative to MSY.

- `B_B0`:

  Time series of depletion.

- `SSB`:

  Time series of spawning stock biomass.

- `SSB_SSBMSY`:

  Time series of spawning stock biomass relative to MSY.

- `SSB_SSB0`:

  Time series of spawning stock depletion.

- `VB`:

  Time series of vulnerable biomass.

- `VB_VBMSY`:

  Time series of vulnerable biomass relative to MSY.

- `VB_VB0`:

  Time series of vulnerable biomass depletion.

- `R`:

  Time series of recruitment.

- `N`:

  Time series of population abundance.

- `N_at_age`:

  Time series of numbers-at-age matrix.

- `Selectivity`:

  Selectivity-at-age matrix.

- `Obs_Catch`:

  Observed catch.

- `Obs_Index`:

  Observed index.

- `Obs_C_at_age`:

  Observed catch-at-age matrix.

- `Catch`:

  Predicted catch.

- `Index`:

  Predicted index.

- `C_at_age`:

  Predicted catch-at-age matrix.

- `Dev`:

  A vector of estimated deviation parameters.

- `Dev_type`:

  A description of the deviation parameters, e.g. "log recruitment
  deviations".

- `NLL`:

  Negative log-likelihood. A vector for the total likelihood, integrated
  across random effects if applicable, components, and penalty term
  (applied when `U > 0.975` in any year).

- `SE_UMSY`:

  Standard error of UMSY estimate.

- `SE_FMSY`:

  Standard error of FMSY estimate.

- `SE_MSY`:

  Standard error of MSY estimate.

- `SE_U_UMSY`:

  Standard error of U/UMSY.

- `SE_F_FMSY`:

  Standard error of F/FMSY.

- `SE_B_BMSY`:

  Standard error of B/BMSY.

- `SE_B_B0`:

  Standard error of B/B0.

- `SE_SSB_SSBMSY`:

  Standard error of SSB/SSBMSY.

- `SE_SSB_SSB0`:

  Standard error of SSB/SSB0.

- `SE_VB_VBMSY`:

  Standard error of VB/VBMSY.

- `SE_VB_VB0`:

  Standard error of VB/VB0.

- `SE_Dev`:

  A vector of standard errors of the deviation parameters.

- `info`:

  A list containing the data and starting values of estimated parameters
  for the assessment.

- `forecast`:

  A list containing components for forecasting:

  - `per_recruit` A data frame of SPR (spawning potential ratio) and YPR
    (yield-per-recruit), calculated for a range of exploitation rate of
    0 - 0.99 or instantaneous F from 0 - 2.5 FMSY.

  - `catch_eq` A function that calculates the catch for the next year
    (after the model terminal year) when an apical F is provided.

- `obj`:

  A list with components returned from
  [`TMB::MakeADFun()`](https://rdrr.io/pkg/TMB/man/MakeADFun.html).

- `opt`:

  A list with components from calling
  [`stats::nlminb()`](https://rdrr.io/r/stats/nlminb.html) to `obj`.

- `SD`:

  A list (class sdreport) with parameter estimates and their standard
  errors, obtained from
  [`TMB::sdreport()`](https://rdrr.io/pkg/TMB/man/sdreport.html).

- `TMB_report`:

  A list of model output reported from the TMB executable, i.e.
  `obj$report()`, and derived quantities (e.g. MSY).

- `dependencies`:

  A character string of data types required for the assessment.

## See also

[plot.Assessment](https://samtool.openmse.com/reference/plot.Assessment.md)
[summary.Assessment](https://samtool.openmse.com/reference/summary.Assessment.md)
[retrospective](https://samtool.openmse.com/reference/retrospective.md)
[profile](https://samtool.openmse.com/reference/profile.md)
[make_MP](https://samtool.openmse.com/reference/make_MP.md)

## Author

Q. Huynh

## Examples

``` r
# \donttest{
output <- DD_TMB(Data = MSEtool::SimulatedData)
class(output)
#> [1] "Assessment"
#> attr(,"package")
#> [1] "SAMtool"
# }
```
