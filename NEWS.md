The latest release of the SAMtool package is available on [CRAN](https://CRAN.R-project.org/package=SAMtool).

## SAMtool 1.2.4
- Update `Shortcut` indexing to align year of assessment with projection. An MP using the `Perfect` assessment and `HCR_MSY` annually will produce F = FMSY in the OM.

## SAMtool 1.2.3
- `RCM` will check age and length comp data for NA's and replaces with zero
- `RCM` reports annual equilibrium unfished reference points using constant stock recruit alpha and beta 
- Unlink TMB executable when unloading package
- Fix punctuation issue in markdown reports that use knitr v1.36.
- Other minor fixes

## SAMtool 1.2.2
- The `make_interim_MP` function is added to generate MPs that adjust the TAC between periodic assessments using the index.
- An additional posfun to `SP` is added to avoid negative biomass situations.
- Fix sampling of recruitment deviations for projections in `RCM` so that the mean is one in normal space. This error was apparent when autocorrelation was very large.
- Assessment and RCModel objects now save the package version as attributes.
- Shortcut function now reports year-specific reference points for the harvest control rule.
- HCR returns TAC = NA when a control point or target point is negative.

## SAMtool 1.2.1
- `HCR_segment` allows for creating control rules with any number of linear segments.
- Assign index beta = 1 for OMs conditioned by RCM.
- The ESS argument is re-instated for RCM to maintain backwards compatibility.
- Fix assignment of length bins in OM by `RCM`.

## SAMtool 1.2.0
### RCM
- A new S4 object, `RCMdata`, is used to send data to the RCM model, i.e., `RCM(OM, RCMdata)`. For now, backwards compatibility should still be maintained when feeding a data list (used prior to v1.2) to fit the model.
- Internal R and TMB code for RCM has been revised, e.g., reducing interchangeability between the terms 'survey' and 'index' to focus on 'index' as much as possible when maintaining backwards compatibility. 
- The `profile` generic is now available for `RCM` models. Steepness, R0, and final depletion can be profiled.
- Fixed plotting bug in `compare_RCM`.
- Priors for index q in `RCM` are now lognormal instead of normal.
- Annual equilibrium and transitional SPR are now calculated and reported.
- Uneven length bin widths are now supported.
- An example is added to the help documentation and online article, using Pacific cod (courtesy of R. Forrest).

### Assessment models
- Priors on M, steepness, R0, and index q (lognormal, see RCM) can be created SCA, DD, and cDD assessment models.
- Likelihood weights to SCA assessments can now be provided for `Catch`, `CAA`, and `CAL` in addition to `Index` in a named list `LWT`. Backwards compatibility remains to provide `LWT` as a vector for index likelihood weights only. 
- A new SCA assessment model that incorporates density-dependent natural mortality (`SCA_DDM`) is added.
- A new SCA assessment model that fits to length composition (`SCA_CAL`) is added.
- Delay-difference assessments (DD_TMB, DD_SS, cDD, cDD_SS) can now fit to mean weight with argument `MW = TRUE`. The functions will look for mean weight data series in `Data@Misc[[x]]$MW`, otherwise will convert length composition `Data@CAL` to weights and calculate annual means.
- Delay-difference assessments (DD_TMB and DD_SS) now use an instantaneous F formulation for the catch equation instead of Pope's approximation. The models should now be more robust for high F situations.
- Dynamic SSB0 is calculated for all assessment models. 
- A variant of the shortcut assessment emulator for closed-loop simulation is available (`Shortcut2`). This function fits an SCA assessment and then characterizes the assessment error relative to the operating model using a vector autoregressive (VAR) model. The functions samples the operating model with error predicted from the VAR model for the projection period. This is a useful function to guide the level of error in the shortcut method.

### Harvest control rules
- Additional OCP types in `HCR_ramp` are available to create harvest control rules based on dynamic B0, and F-based rules (F/FMSY, F/F01, F/F-SPR).
- Added a shortcut function for a fixed escapement harvest control rule (`HCR_escapement`).

## SAMtool 1.1.2
- Minor fix to vignette to fix MSEtool reverse dependency issue.

## SAMtool 1.1.1
- Fix year range for depletion calculation from time-varying SSB0 in RCM.

## SAMtool 1.1.0
- Edits to fix valgrind and sanitizer issues in TMB code.
- Likelihood gradients, the derivatives of the likelihood of each annual data point with respect to model parameters, are plotted in the RCM markdown report. For this purpose, the annual age or length composition is considered to be a single piece of data. This diagnostic could be informative on how informative the data are to model parameters, with more influential data points having larger gradients.
- When conditioned on effort, the `RCM` will now incorporate catches into the likelihood as a default. This allows the model to estimate F and R0 when conditioned on effort and there is patchy catch data.

## SAMtool 1.0.0
- Development of the assessment models and OM conditioning model in SAMtool 1.0.0 continues from MSEtool 2.0.1. `multiMSE` remains in MSEtool.
- The age structure of the SCA models (`SCA`, `SCA_Pope`, `SSS`) start at age 0 following the change in the MSEtool OM.
- An additional SCA model (`SCA_RWM`) can be used to estimate time-varying M (constant with age) as a random walk. Fix the random walk SD to a low value to effectively estimate a time-constant M (see help page).
- Warnings during the fit of the assessment models (through `nlminb`) are turned off. Convergence status and issues can be checked in the `conv` slot of the output Assessment object. In closed-loop simulation, the `diagnostic` function can be used to track the behavior of model-based MPs. By default, pre-packaged model-based MPs and MPs made from the `make_MP` function are designed to report convergence info (stored in `MSE@PPD`). 
- Assessment functions now calculate and report spawning potential ratio and yield per recruit in the forecast slot of the S4 object. Also in this slot is a catch equation function calculates the TAC for a given F. 
- HCR nomenclature has changed. Operational control points (OCPs) are used instead of reference points (to help distinguish between reference points in the estimation model vs. the operating model. Various types of F can now be used in the HCR, including F0.1, Fmax, and FSPR, in addition to FMSY for the TAC calculation.
- A `Shortcut` assess function samples the OM with error and autocorrelation for HCRs as an emulator of a stock assessment in closed-loop simulation. The `Perfect` function samples the OM without error.
- All assessment models now accommodate multiple indices in model fitting, specify in the `AddInd` argument of functions which index slots in the Data object will be used among Data@Ind, Data@SpInd, Data@VInd, and Data@AddInd. Within series weighting is applied by using the corresponding CV slot, i.e., Data@CV_Ind for Data@Ind, etc. Among series weighting can also be tuned using likelihood weights with `LWT` argument. For SCA and VPA models, the selectivity is fixed in the model using Data@AddIndV for indices in Data@AddInd. 

### RCM
- The function for the OM conditioning model is now re-named to `RCM` (Rapid Conditioning Model). 
- The age structure of the model now starts at age 0, following the change in the MSEtool OM. The dimension associated with age in matrices and arrays need to be of length `maxage + 1` which corresponds to ages 0 to maxage.
- When estimating fleet F in the model (`condition = "catch"`), the likelihood for the catch can now have a user-defined standard deviation indicated in `data$C_sd` (year and fleet specific, the previous default was 0.01 was built-in for all catches).
- For generating length comps, variability in length-at-age can now be age-specific. Specify the length-at-age standard deviation in `OM@cpars$LatASD`.
- Priors for log_R0 (normal distribution), steepness (beta distribution for Beverton-Holt, normal for Ricker), log_M (age and time constant, normal), and survey q (normal) can now be specified.
