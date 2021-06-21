The latest release of the SAMtool package is available on [CRAN](https://CRAN.R-project.org/package=SAMtool).

## SAMtool 1.2.0
### RCM
- A new S4 object, `RCMdata`, used to send data to the RCM model, i.e., `RCM(OM, RCMdata)`. Internal R and TMB code for RCM has been revised, e.g., reducing interchangeability between the terms 'survey' and 'index' to focus on 'index' as much as possible when maintaining backwards compatibility. For now, backwards compatibility should still be maintained when feeding a data list (used prior to v1.2) to fit the model.
- Create a profiling function for `RCM` models. Steepness, R0, and final depletion can be profiled.
- Fix plotting bug in `compare_RCM`.
- Priors for index q in `RCM` is now lognormal instead of normal.
- Annual equilibrium and transitional SPR are now calculated and reported.
- Uneven length bin widths now supported.
- An example is now added to the help documentation and online article, using Pacific cod (courtesy of R. Forrest).

### Assessment models
- Priors on M, steepness, R0, and index q (lognormal, see RCM) can be created SCA, DD, and cDD assessment models.
- Likelihood weights to SCA assessments can now be provided for `Catch`, `CAA`, and `CAL` in addition to `Index` in a named list `LWT`.
Backwards compatibility remains to provide `LWT` as a vector for indices weights only. 
- A new SCA assessment model that incorporates density-dependent natural mortality (`SCA_DDM`) is added.
- A new SCA assessment model that fits to length composition (`SCA_CAL`) is added.
- Delay-difference assessments (DD_TMB, DD_SS, cDD, cDD_SS) can now fit to mean weight with argument `MW = TRUE`. The functions will look for mean weight data series in `Data@Misc[[x]]$MW`, otherwise will convert length composition `Data@CAL` to weights and calculate annual means.
- Delay-difference assessments (DD_TMB and DD_SS) use an instantaneous F formulation for the catch equation instead of Pope's approximation. The models should now be more robust for high F situations.
- Dynamic SSB0 is calculated for all assessment models. 

### Harvest control rules
- Additional OCP types in `HCR_ramp` to create harvest control rules based on dynamic B0, and F-based rules (F/FMSY, F/F01, F/F-SPR).
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
