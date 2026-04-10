# Assessment emulator as a shortcut to model fitting in closed-loop simulation

Functions (class Assessment) that emulate a stock assessment by sampling
the operating model biomass, abundance, and fishing mortality (with
observation error, autocorrelation, and bias) instead of fitting a
model. This output can then be passed onto a harvest control rule (HCR
function). `Shortcut` is the base function that samples the OM with an
error distribution. `Shortcut2`, the more preferable option, fits
[SCA](https://samtool.openmse.com/reference/SCA.md) in the last
historical year of the operating model, estimates the error parameters
using a vector autoregressive model of the residuals, and then generates
model "estimates" using
[predict.varest](https://rdrr.io/pkg/vars/man/predict.html). `Perfect`
assumes no error in the assessment model and is useful for comparing the
behavior of different harvest control rules. To utilize the shortcut
method in closed-loop simulation, use
[make_MP](https://samtool.openmse.com/reference/make_MP.md) with these
functions as the Assessment model. **N.B. the functions do not work
with** `runMSE(parallel = TRUE)` for MSEtool v3.4.0 and earlier.

## Usage

``` r
Shortcut(
  x = 1,
  Data,
  method = c("B", "N", "RF"),
  B_err = c(0.3, 0.7, 1),
  N_err = c(0.3, 0.7, 1),
  R_err = c(0.3, 0.7, 1),
  F_err = c(0.3, 0.7, 1),
  VAR_model,
  ...
)

Shortcut2(
  x,
  Data,
  method = "N",
  SCA_args = list(),
  VAR_args = list(type = "none"),
  ...
)

Perfect(x, Data, ...)
```

## Arguments

- x:

  An index for the objects in `Data` when running in
  [runMSE](https://msetool.openmse.com/reference/runMSE.html).
  Otherwise, equals to 1 When running an assessment interactively.

- Data:

  An object of class Data.

- method:

  Indicates where the error in the OM is located. For "B", OM biomass is
  directly sampled with error. For "N", OM abundance-at-age is sampled
  and biomass subsequently calculated. For "RF", recruitment and F are
  sampled to calculate abundance and biomass. There is no error in
  biological parameters for "N" and "RF". By default, "B" is used for
  `Shortcut` and "N" for `Shortcut2`.

- B_err:

  If `method = "B"`, a vector of length three that specifies the
  standard deviation (in logspace), autocorrelation, and bias (1 =
  unbiased) for biomass.

- N_err:

  Same as B_err, but for abundance when `method = "N"`.

- R_err:

  Same as B_err, but for recruitment when `method = "RF"`.

- F_err:

  Same as B_err. Always used regardless of `method` to report F and
  selectivity for HCR.

- VAR_model:

  An object returned by [VAR](https://rdrr.io/pkg/vars/man/VAR.html) to
  generate emulated assessment error. Used by `Shortcut2`.

- ...:

  Other arguments (not currently used).

- SCA_args:

  Additional arguments to pass to
  [SCA](https://samtool.openmse.com/reference/SCA.md). Currently,
  arguments `SR` and `vulnerability` are obtained from the operating
  model.

- VAR_args:

  Additional arguments to pass to
  [VAR](https://rdrr.io/pkg/vars/man/VAR.html). By default, argument
  `type = "none"` (stationary time series with mean zero is assumed).

## Value

An object of class
[Assessment](https://samtool.openmse.com/reference/Assessment-class.md).

## Details

Currently there is no error in FMSY (frequently the target F in the
HCR).

See Wiedenmann et al. (2015) for guidance on the magnitude of error for
the shortcut emulator.

## References

Wiedenmann, J., Wilberg, M.J., Sylvia, A., and Miller, T.J. 2015.
Autocorrelated error in stock assessment estimates: Implications for
management strategy evaluation. Fisheries Research 172: 325-334.

## Author

Q. Huynh

## Examples

``` r
Shortcut_4010 <- make_MP(Shortcut, HCR40_10) 
Shortcut_Nerr <- make_MP(Shortcut, HCR40_10, method = "N", N_err = c(0.1, 0.1, 1)) # Highly precise!

# Fits SCA first and then emulate it in the projection period 
Shortcut2_4010 <- make_MP(Shortcut2, HCR40_10) 

# \donttest{
# Compare the shortcut method vs. fitting an SCA model with a 40-10 control rule
MSE <- runMSE(testOM, MPs = c("Shortcut_4010", "SCA_4010"))
#> → Loading operating model
#> → Calculating MSY reference points for each year
#> → Optimizing for user-specified depletion in last historical year
#> → Calculating historical stock and fishing dynamics
#> → Calculating per-recruit reference points
#> → Calculating B-low reference points
#> → Calculating reference yield - best fixed F strategy
#> → Simulating observed data
#> ✔ Running forward projections
#> 
#> → 1/2 Running MSE for: "Shortcut_4010" 
#>   |==                                                | 2 % ~01s            |===                                               | 4 % ~01s            |====                                              | 6 % ~01s            |=====                                             | 8 % ~01s            |======                                            | 10% ~01s            |=======                                           | 12% ~01s            |========                                          | 14% ~01s            |=========                                         | 16% ~01s            |==========                                        | 18% ~01s            |===========                                       | 20% ~01s            |============                                      | 22% ~01s            |=============                                     | 24% ~01s            |==============                                    | 27% ~01s            |===============                                   | 29% ~01s            |================                                  | 31% ~01s            |=================                                 | 33% ~01s            |==================                                | 35% ~01s            |===================                               | 37% ~01s            |====================                              | 39% ~01s            |=====================                             | 41% ~01s            |======================                            | 43% ~01s            |=======================                           | 45% ~01s            |========================                          | 47% ~01s            |=========================                         | 49% ~01s            |==========================                        | 51% ~01s            |===========================                       | 53% ~01s            |============================                      | 55% ~01s            |=============================                     | 57% ~01s            |==============================                    | 59% ~01s            |===============================                   | 61% ~01s            |================================                  | 63% ~00s            |=================================                 | 65% ~00s            |==================================                | 67% ~00s            |===================================               | 69% ~00s            |====================================              | 71% ~00s            |=====================================             | 73% ~00s            |======================================            | 76% ~00s            |=======================================           | 78% ~00s            |========================================          | 80% ~00s            |=========================================         | 82% ~00s            |==========================================        | 84% ~00s            |===========================================       | 86% ~00s            |============================================      | 88% ~00s            |=============================================     | 90% ~00s            |==============================================    | 92% ~00s            |===============================================   | 94% ~00s            |================================================  | 96% ~00s            |================================================= | 98% ~00s            |==================================================| 100% elapsed=01s  
#> 
#> 
#> → 2/2 Running MSE for: "SCA_4010" 
#>   |==                                                | 2 % ~54s            |===                                               | 4 % ~27s            |====                                              | 6 % ~18s            |=====                                             | 8 % ~13s            |======                                            | 10% ~22s            |=======                                           | 12% ~18s            |========                                          | 14% ~15s            |=========                                         | 16% ~13s            |==========                                        | 18% ~17s            |===========                                       | 20% ~15s            |============                                      | 22% ~14s            |=============                                     | 24% ~12s            |==============                                    | 27% ~15s            |===============                                   | 29% ~14s            |================                                  | 31% ~12s            |=================                                 | 33% ~11s            |==================                                | 35% ~14s            |===================                               | 37% ~13s            |====================                              | 39% ~12s            |=====================                             | 41% ~11s            |======================                            | 43% ~12s            |=======================                           | 45% ~11s            |========================                          | 47% ~11s            |=========================                         | 49% ~10s            |==========================                        | 51% ~11s            |===========================                       | 53% ~10s            |============================                      | 55% ~09s            |=============================                     | 57% ~09s            |==============================                    | 59% ~09s            |===============================                   | 61% ~09s            |================================                  | 63% ~08s            |=================================                 | 65% ~07s            |==================================                | 67% ~08s            |===================================               | 69% ~07s            |====================================              | 71% ~07s            |=====================================             | 73% ~06s            |======================================            | 76% ~06s            |=======================================           | 78% ~06s            |========================================          | 80% ~05s            |=========================================         | 82% ~04s            |==========================================        | 84% ~05s            |===========================================       | 86% ~04s            |============================================      | 88% ~03s            |=============================================     | 90% ~03s            |==============================================    | 92% ~02s            |===============================================   | 94% ~02s            |================================================  | 96% ~01s            |================================================= | 98% ~01s            |==================================================| 100% elapsed=31s  
#> 
# }

# Compare the performance of three HCRs
Perfect_4010 <- make_MP(Perfect, HCR40_10)
Perfect_6020 <- make_MP(Perfect, HCR60_20)
Perfect_8040MSY <- make_MP(Perfect, HCR_ramp, OCP_type = "SSB_SSBMSY", TOCP = 0.8, LOCP = 0.4)

# \donttest{
MSE <- runMSE(testOM, MPs = c("Perfect_4010", "Perfect_6020", "Perfect_8040MSY"))
#> → Loading operating model
#> → Calculating MSY reference points for each year
#> → Optimizing for user-specified depletion in last historical year
#> → Calculating historical stock and fishing dynamics
#> → Calculating per-recruit reference points
#> → Calculating B-low reference points
#> → Calculating reference yield - best fixed F strategy
#> → Simulating observed data
#> ✔ Running forward projections
#> 
#> → 1/3 Running MSE for: "Perfect_4010" 
#>   |==                                                | 2 % ~01s            |===                                               | 4 % ~01s            |====                                              | 6 % ~00s            |=====                                             | 8 % ~00s            |======                                            | 10% ~01s            |=======                                           | 12% ~01s            |========                                          | 14% ~01s            |=========                                         | 16% ~00s            |==========                                        | 18% ~01s            |===========                                       | 20% ~01s            |============                                      | 22% ~01s            |=============                                     | 24% ~01s            |==============                                    | 27% ~01s            |===============                                   | 29% ~01s            |================                                  | 31% ~01s            |=================                                 | 33% ~00s            |==================                                | 35% ~01s            |===================                               | 37% ~01s            |====================                              | 39% ~00s            |=====================                             | 41% ~00s            |======================                            | 43% ~00s            |=======================                           | 45% ~00s            |========================                          | 47% ~00s            |=========================                         | 49% ~00s            |==========================                        | 51% ~00s            |===========================                       | 53% ~00s            |============================                      | 55% ~00s            |=============================                     | 57% ~00s            |==============================                    | 59% ~00s            |===============================                   | 61% ~00s            |================================                  | 63% ~00s            |=================================                 | 65% ~00s            |==================================                | 67% ~00s            |===================================               | 69% ~00s            |====================================              | 71% ~00s            |=====================================             | 73% ~00s            |======================================            | 76% ~00s            |=======================================           | 78% ~00s            |========================================          | 80% ~00s            |=========================================         | 82% ~00s            |==========================================        | 84% ~00s            |===========================================       | 86% ~00s            |============================================      | 88% ~00s            |=============================================     | 90% ~00s            |==============================================    | 92% ~00s            |===============================================   | 94% ~00s            |================================================  | 96% ~00s            |================================================= | 98% ~00s            |==================================================| 100% elapsed=01s  
#> 
#> 
#> → 2/3 Running MSE for: "Perfect_6020" 
#>   |==                                                | 2 % ~01s            |===                                               | 4 % ~01s            |====                                              | 6 % ~00s            |=====                                             | 8 % ~00s            |======                                            | 10% ~01s            |=======                                           | 12% ~01s            |========                                          | 14% ~01s            |=========                                         | 16% ~00s            |==========                                        | 18% ~01s            |===========                                       | 20% ~01s            |============                                      | 22% ~01s            |=============                                     | 24% ~00s            |==============                                    | 27% ~01s            |===============                                   | 29% ~01s            |================                                  | 31% ~01s            |=================                                 | 33% ~00s            |==================                                | 35% ~01s            |===================                               | 37% ~00s            |====================                              | 39% ~00s            |=====================                             | 41% ~00s            |======================                            | 43% ~00s            |=======================                           | 45% ~00s            |========================                          | 47% ~00s            |=========================                         | 49% ~00s            |==========================                        | 51% ~00s            |===========================                       | 53% ~00s            |============================                      | 55% ~00s            |=============================                     | 57% ~00s            |==============================                    | 59% ~00s            |===============================                   | 61% ~00s            |================================                  | 63% ~00s            |=================================                 | 65% ~00s            |==================================                | 67% ~00s            |===================================               | 69% ~00s            |====================================              | 71% ~00s            |=====================================             | 73% ~00s            |======================================            | 76% ~00s            |=======================================           | 78% ~00s            |========================================          | 80% ~00s            |=========================================         | 82% ~00s            |==========================================        | 84% ~00s            |===========================================       | 86% ~00s            |============================================      | 88% ~00s            |=============================================     | 90% ~00s            |==============================================    | 92% ~00s            |===============================================   | 94% ~00s            |================================================  | 96% ~00s            |================================================= | 98% ~00s            |==================================================| 100% elapsed=01s  
#> 
#> 
#> → 3/3 Running MSE for: "Perfect_8040MSY" 
#>   |==                                                | 2 % ~01s            |===                                               | 4 % ~01s            |====                                              | 6 % ~00s            |=====                                             | 8 % ~00s            |======                                            | 10% ~01s            |=======                                           | 12% ~01s            |========                                          | 14% ~01s            |=========                                         | 16% ~01s            |==========                                        | 18% ~01s            |===========                                       | 20% ~01s            |============                                      | 22% ~01s            |=============                                     | 24% ~01s            |==============                                    | 27% ~01s            |===============                                   | 29% ~01s            |================                                  | 31% ~01s            |=================                                 | 33% ~00s            |==================                                | 35% ~01s            |===================                               | 37% ~00s            |====================                              | 39% ~00s            |=====================                             | 41% ~00s            |======================                            | 43% ~00s            |=======================                           | 45% ~00s            |========================                          | 47% ~00s            |=========================                         | 49% ~00s            |==========================                        | 51% ~00s            |===========================                       | 53% ~00s            |============================                      | 55% ~00s            |=============================                     | 57% ~00s            |==============================                    | 59% ~00s            |===============================                   | 61% ~00s            |================================                  | 63% ~00s            |=================================                 | 65% ~00s            |==================================                | 67% ~00s            |===================================               | 69% ~00s            |====================================              | 71% ~00s            |=====================================             | 73% ~00s            |======================================            | 76% ~00s            |=======================================           | 78% ~00s            |========================================          | 80% ~00s            |=========================================         | 82% ~00s            |==========================================        | 84% ~00s            |===========================================       | 86% ~00s            |============================================      | 88% ~00s            |=============================================     | 90% ~00s            |==============================================    | 92% ~00s            |===============================================   | 94% ~00s            |================================================  | 96% ~00s            |================================================= | 98% ~00s            |==================================================| 100% elapsed=01s  
#> 
# }
```
