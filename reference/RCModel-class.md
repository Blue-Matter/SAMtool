# Class-`RCModel`

An S4 class for the output from
[RCM](https://samtool.openmse.com/reference/RCM.md).

## Slots

- `OM`:

  An updated operating model, class
  [MSEtool::OM](https://msetool.openmse.com/reference/OM-class.html).

- `SSB`:

  A matrix of estimated spawning biomass with `OM@nsim` rows and
  `OM@nyears+1` columns.

- `NAA`:

  An array for the predicted numbers at age with dimension `OM@nsim`,
  `OM@nyears+1`, and `OM@maxage+1`.

- `CAA`:

  An array for the predicted catch at age with dimension `OM@nsim`,
  `OM@nyears`, `OM@maxage`, and nfleet.

- `CAL`:

  An array for the predicted catch at length with dimension `OM@nsim`,
  `OM@nyears`, length bins, and nfleet.

- `conv`:

  A logical vector of length `OM@nsim` indicating convergence of the RCM
  in the i-th simulation.

- `report`:

  A list of length `OM@nsim` with more output from the fitted RCM.
  Within each simulation, a named list containing items of interest
  include:

  - B - total biomass - vector of length nyears+1

  - EPR0 - annual unfished spawners per recruit - vector of length
    nyears

  - ageM - age of 50% maturity - integer

  - EPR0_SR - unfished spawners per recruit for the stock-recruit
    relationship (mean EPR0 over the first `ageM` years) - numeric

  - R0 - unfished recruitment for the stock-recruit relationship -
    numeric

  - h - steepness for the stock-recruit relationship - numeric

  - Arec - stock-recruit alpha parameter - numeric

  - Brec - stock-recruit beta parameter - numeric

  - E0_SR - unfished spawning biomass for the stock-recruit relationship
    (product of EPR0_SR and R0) - numeric

  - CR_SR - compensation ratio, the product of Arec and EPR0_SR -
    numeric

  - E0 - annual unfished spawning biomass (intersection of stock-recruit
    relationship and unfished spawners per recruit) - vector of length
    nyears

  - R0_annual - annual unfished recruitment (annual ratio of E0 and
    EPR0) - vector of length nyears

  - h_annual - annual steepness (calculated from EPR0 and Arec) - vector
    of length nyears

  - CR - annual compensation ratio, the product of alpha and annual
    unfished spawners per recruit (EPR0) - vector of length nyears

  - R - recruitment - vector of length nyears+1

  - R_early - recruitment for the cohorts in first year of the model -
    vector n_age-1 (where n_age = maxage + 1)

  - VB - vulnerable biomass - matrix of nyears x nfleet

  - N - abundance at age - matrix of nyears+1 x n_age

  - F - apical fishing mortality - matrix of nyears x nfleet

  - F_at_age - fishing mortality at age - matrix of nyears x n_age

  - F_equilibrium - equilibrium fishing mortality prior to first year -
    vector of length nfleet

  - M - natural mortality - matrix of nyears x n_age

  - Z - total mortality - matrix of nyears x n_age

  - q - index catchability - vector of length nsurvey

  - ivul - index selectivity at age - array of dim nyears+1, n_age,
    nsurvey

  - ivul_len - corresponding index selectivity at length - matrix of
    nbins x nsurvey

  - Ipred - predicted index values - matrix of nyears x nsurvey

  - IAApred - predicted index catch at age - array of dim nyears, n_age,
    nsurvey

  - vul - fleet selectivity at age - array of dim nyears+1, n_age,
    nfleet (or nsel_block)

  - vul_len - corresponding fleet selectivity at length - matrix of
    nbins x nfleet (or nsel_block)

  - IALpred - predicted index catch at length - array of dim nyears,
    nbins, nsurvey

  - MLpred - predicted mean length - matrix of nyears x nfleet

  - MWpred - predicted mean weight - matrix of nyears x nfleet

  - CAApred - predicted catch at age - array of nyears, n_age, nfleet

  - CALpred - predicted catch at length - array of nyears, nbins, nfleet

  - Cpred - predicted catch in weight - matrix of nyears x nfleet

  - CN - predicted catch in numbers - matrix of nyears x nfleet

  - dynamic_SSB0 - the dynamic unfished spawning biomass calcaluated by
    projecting the historical model with zero catches - vector of length
    nyears+1

  - SPR_eq - equilibrium spawning potential ratio calculated from annual
    F-at-age - vector of length nyears

  - SPR_dyn - dynamic (transitional) spawning potential ratio calculated
    from cumulative survival of cohorts - vector of length nyears

  - nll - total objective function of the model - numeric

  - nll_fleet - objective function values for each annual data point(s)
    from fleets - array of nyears x nfleet x 5 (for Catch, equilibrium
    catch, CAA, CAL, and mean size)

  - nll_index - objective function values for each annual data point(s)
    in the index - array of nyears x nsurvey x 3 (for Index, IAA, and
    IAL)

  - prior - penalty value added to the objective function from priors -
    numeric

  - penalty - additional penalty values added to the objective function
    due to high F - numeric

  - conv - whether the model converged (whether a positive-definite
    Hessian was obtained) - logical

- `mean_fit`:

  A list of output from fit to mean values of life history parameters in
  the operating model. The named list consists of:

  - obj - a list with components returned from
    [`TMB::MakeADFun()`](https://rdrr.io/pkg/TMB/man/MakeADFun.html).

  - opt - a list with components from calling
    [`stats::nlminb()`](https://rdrr.io/r/stats/nlminb.html) to `obj`.

  - SD - a list (class sdreport) with parameter estimates and their
    standard errors, obtained from
    [`TMB::sdreport()`](https://rdrr.io/pkg/TMB/man/sdreport.html).

  - report - a list of model output reported from the TMB executable,
    i.e. `obj$report()`. See Misc.

- `data`:

  A [RCMdata](https://samtool.openmse.com/reference/RCMdata-class.md)
  object containing data inputs for the RCM.

- `config`:

  A list describing configuration of the RCM:

  - drop_sim - a vector of simulations that were dropped for the output

- `Misc`:

  Slot for miscellaneous information for the user. Currently unused.

## See also

[plot.RCModel](https://samtool.openmse.com/reference/plot.RCModel.md)
[RCM](https://samtool.openmse.com/reference/RCM.md)

## Author

Q. Huynh
