# Fixed escapement harvest control rule

A simple control rule that allows fishing when the operational control
point (OCP) is above some threshold. By default, this function sets the
TAC at F = 100% FMSY when spawning depletion \> 0.1.

## Usage

``` r
HCR_escapement(
  Assessment,
  reps = 1,
  OCP_type = "SSB_SSB0",
  OCP_threshold = 0.2,
  Ftarget_type = "FMSY",
  relF_max = 1,
  ...
)
```

## Arguments

- Assessment:

  An object of class
  [Assessment](https://samtool.openmse.com/reference/Assessment-class.md)
  with estimates of FMSY or UMSY and vulnerable biomass in terminal
  year.

- reps:

  The number of stochastic samples of the TAC recommendation.

- OCP_type:

  The type of operational control points (OCPs) for the harvest control
  rule used to determine whether there is fishing. By default, use
  (`"SSB_SSB0"` for spawning depletion. Other biomass OCPs include
  `"SSB_SSBMSY"` for spawning biomass relative to MSY and `"SSB_dSSB0"`,
  for dynamic depletion (dynamic SSB0 is the historical reconstructed
  biomass with F = 0). For F-based OCPs, the terminal year fishing
  mortality relative F01 or Fmax (using yield-per-recruit) or F-SPR%
  (see `SPR_OCP` argument) can be used.

- OCP_threshold:

  The value of the OCP above which fishing can occur.

- Ftarget_type:

  The type of F used for the target fishing mortality rate.

- relF_max:

  The relative value of Ftarget if `OCP > OCP_treshold`.

- ...:

  Miscellaneous arguments.

## Value

An object of class
[MSEtool::Rec](https://msetool.openmse.com/reference/Rec-class.html)
with the TAC recommendation.

## Details

The catch advice is calculated using the catch equation of the
corresponding assessment. See `Assessment@forecast$catch_eq`, a function
that returns the catch advice for a specified `Ftarget`.

## References

Deroba, J.J. and Bence, J.R. 2008. A review of harvest policies:
Understanding relative performance of control rules. Fisheries Research
94:210-223.

## See also

[make_MP](https://samtool.openmse.com/reference/make_MP.md)
[HCR_ramp](https://samtool.openmse.com/reference/HCR_ramp.md)

## Author

Q. Huynh

## Examples

``` r
# create an MP to run in closed-loop MSE (fishes at FMSY when B/B0 > 0.2)
SP_escapement <- make_MP(SP, HCR_escapement)

# The MP which fishes at 75% of FMSY
SP_escapement75 <- make_MP(SP, HCR_escapement, relF_max = 0.75)

# The MP which fishes at FMSY when BMSY > 0.5
SP_BMSY_escapement <- make_MP(SP, HCR_escapement, OCP_type = "SSB_SSBMSY", 
                              OCP_threshold = 0.5, relF_max = 1)

# \donttest{
myOM <- MSEtool::runMSE(MSEtool::testOM, MPs = c("FMSYref", "SP_escapement", "SP_BMSY_escapement"))
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
#> → 1/3 Running MSE for: "FMSYref" 
#>   |==                                                | 2 % ~00s            |===                                               | 4 % ~00s            |====                                              | 6 % ~00s            |=====                                             | 8 % ~00s            |======                                            | 10% ~01s            |=======                                           | 12% ~00s            |========                                          | 14% ~00s            |=========                                         | 16% ~00s            |==========                                        | 18% ~00s            |===========                                       | 20% ~00s            |============                                      | 22% ~00s            |=============                                     | 24% ~00s            |==============                                    | 27% ~00s            |===============                                   | 29% ~00s            |================                                  | 31% ~00s            |=================                                 | 33% ~00s            |==================                                | 35% ~00s            |===================                               | 37% ~00s            |====================                              | 39% ~00s            |=====================                             | 41% ~00s            |======================                            | 43% ~00s            |=======================                           | 45% ~00s            |========================                          | 47% ~00s            |=========================                         | 49% ~00s            |==========================                        | 51% ~00s            |===========================                       | 53% ~00s            |============================                      | 55% ~00s            |=============================                     | 57% ~00s            |==============================                    | 59% ~00s            |===============================                   | 61% ~00s            |================================                  | 63% ~00s            |=================================                 | 65% ~00s            |==================================                | 67% ~00s            |===================================               | 69% ~00s            |====================================              | 71% ~00s            |=====================================             | 73% ~00s            |======================================            | 76% ~00s            |=======================================           | 78% ~00s            |========================================          | 80% ~00s            |=========================================         | 82% ~00s            |==========================================        | 84% ~00s            |===========================================       | 86% ~00s            |============================================      | 88% ~00s            |=============================================     | 90% ~00s            |==============================================    | 92% ~00s            |===============================================   | 94% ~00s            |================================================  | 96% ~00s            |================================================= | 98% ~00s            |==================================================| 100% elapsed=01s  
#> 
#> 
#> → 2/3 Running MSE for: "SP_escapement" 
#>   |==                                                | 2 % ~03s            |===                                               | 4 % ~02s            |====                                              | 6 % ~01s            |=====                                             | 8 % ~01s            |======                                            | 10% ~02s            |=======                                           | 12% ~01s            |========                                          | 14% ~01s            |=========                                         | 16% ~01s            |==========                                        | 18% ~01s            |===========                                       | 20% ~01s            |============                                      | 22% ~01s            |=============                                     | 24% ~01s            |==============                                    | 27% ~01s            |===============                                   | 29% ~01s            |================                                  | 31% ~01s            |=================                                 | 33% ~01s            |==================                                | 35% ~01s            |===================                               | 37% ~01s            |====================                              | 39% ~01s            |=====================                             | 41% ~01s            |======================                            | 43% ~01s            |=======================                           | 45% ~01s            |========================                          | 47% ~01s            |=========================                         | 49% ~01s            |==========================                        | 51% ~01s            |===========================                       | 53% ~01s            |============================                      | 55% ~01s            |=============================                     | 57% ~01s            |==============================                    | 59% ~01s            |===============================                   | 61% ~01s            |================================                  | 63% ~01s            |=================================                 | 65% ~00s            |==================================                | 67% ~00s            |===================================               | 69% ~00s            |====================================              | 71% ~00s            |=====================================             | 73% ~00s            |======================================            | 76% ~00s            |=======================================           | 78% ~00s            |========================================          | 80% ~00s            |=========================================         | 82% ~00s            |==========================================        | 84% ~00s            |===========================================       | 86% ~00s            |============================================      | 88% ~00s            |=============================================     | 90% ~00s            |==============================================    | 92% ~00s            |===============================================   | 94% ~00s            |================================================  | 96% ~00s            |================================================= | 98% ~00s            |==================================================| 100% elapsed=02s  
#> 
#> 
#> → 3/3 Running MSE for: "SP_BMSY_escapement" 
#>   |==                                                | 2 % ~03s            |===                                               | 4 % ~02s            |====                                              | 6 % ~01s            |=====                                             | 8 % ~01s            |======                                            | 10% ~02s            |=======                                           | 12% ~01s            |========                                          | 14% ~01s            |=========                                         | 16% ~01s            |==========                                        | 18% ~01s            |===========                                       | 20% ~01s            |============                                      | 22% ~01s            |=============                                     | 24% ~01s            |==============                                    | 27% ~01s            |===============                                   | 29% ~01s            |================                                  | 31% ~01s            |=================                                 | 33% ~01s            |==================                                | 35% ~01s            |===================                               | 37% ~01s            |====================                              | 39% ~01s            |=====================                             | 41% ~01s            |======================                            | 43% ~01s            |=======================                           | 45% ~01s            |========================                          | 47% ~01s            |=========================                         | 49% ~01s            |==========================                        | 51% ~01s            |===========================                       | 53% ~01s            |============================                      | 55% ~01s            |=============================                     | 57% ~01s            |==============================                    | 59% ~01s            |===============================                   | 61% ~01s            |================================                  | 63% ~01s            |=================================                 | 65% ~00s            |==================================                | 67% ~00s            |===================================               | 69% ~00s            |====================================              | 71% ~00s            |=====================================             | 73% ~00s            |======================================            | 76% ~00s            |=======================================           | 78% ~00s            |========================================          | 80% ~00s            |=========================================         | 82% ~00s            |==========================================        | 84% ~00s            |===========================================       | 86% ~00s            |============================================      | 88% ~00s            |=============================================     | 90% ~00s            |==============================================    | 92% ~00s            |===============================================   | 94% ~00s            |================================================  | 96% ~00s            |================================================= | 98% ~00s            |==================================================| 100% elapsed=02s  
#> 
# }
```
