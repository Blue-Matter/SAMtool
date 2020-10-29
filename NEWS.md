The latest release of the SAMtool package is available on [CRAN](https://CRAN.R-project.org/package=SAMtool).

## SAMtool 1.0.0
Development of SAMtool continues from SAMtool 2.0.1.

### RCM
- The function for the OM conditioning model is now re-named to `RCM` (Rapid Conditioning Model). 
- The age structure of the model now starts at age 0, following the change in MSEtool. The dimension associated with age in matrices and arrays need to be of length `maxage + 1` which corresponds to ages 0 to maxage.
- When estimating fleet F in the model (`condition = "catch"`), the likelihood for the catch can now have a user-defined standard deviation indicated in `data$C_sd` (year and fleet specific, the previous default was 0.01 was built-in for all catches).
- Variability in length-at-age can vary by age. Specify in `OM@cpars$LatASD`.

