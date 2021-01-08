The latest release of the SAMtool package is available on [CRAN](https://CRAN.R-project.org/package=SAMtool).

## SAMtool 1.0.0
- Development of the assessment models and OM conditioning model in SAMtool 1.0.0 continues from MSEtool 2.0.1. `multiMSE` remains in MSEtool.
- The age structure of the SCA models (`SCA`, `SCA_Pope`, `SSS`) start at age 0 following the change in the MSEtool OM.
- An additional SCA model (`SCA_RWM`) can be used to estimate time-varying M (constant with age) as a random walk. Fix the random walk SD to a low value to effectively estimate a time-constant M (see help page).
- Warnings during the fit of the assessment models (through `nlminb`) are turned off. Convergence status and issues can be checked in the `conv` slot of the output Assessment object. In closed-loop simulation, the `diagnostic` function can be used to track the behavior of model-based MPs. By default, pre-packaged model-based MPs and MPs made from the `make_MP` function are designed to report convergence info (stored in `MSE@PPD`). 
- Assessment functions now calculate and report spawning potential ratio and yield per recruit in the forecast slot of the S4 object. Also in this slot is a catch equation function calculates the TAC for a given F. 
- HCR nomenclature has changed. Operational control points (OCPs) are used instead of reference points (to help distinguish between reference points in the estimation model vs. the operating model. Various types of F can now be used in the HCR, including F0.1, Fmax, and FSPR, in addition to FMSY for the TAC calculation.
- A `Shortcut` assess function samples the OM with error and autocorrelation for HCRs as an emulator of a stock assessment in closed-loop simulation. The `Perfect` function samples the OM without error.

### RCM
- The function for the OM conditioning model is now re-named to `RCM` (Rapid Conditioning Model). 
- The age structure of the model now starts at age 0, following the change in the MSEtool OM. The dimension associated with age in matrices and arrays need to be of length `maxage + 1` which corresponds to ages 0 to maxage.
- When estimating fleet F in the model (`condition = "catch"`), the likelihood for the catch can now have a user-defined standard deviation indicated in `data$C_sd` (year and fleet specific, the previous default was 0.01 was built-in for all catches).
- For generating length comps, variability in length-at-age can now be age-specific. Specify the length-at-age standard deviation in `OM@cpars$LatASD`.
- Priors for log_R0 (normal distribution), steepness (beta distribution for Beverton-Holt, normal for Ricker), log_M (age and time constant, normal), and survey q (normal) can now be specified.
