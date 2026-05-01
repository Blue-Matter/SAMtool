# Model-based management procedures

A suite of model-based management procedures (MPs) included in the
package. Additional MPs, with specific model configurations (e.g.,
stock-recruit function or fixing certain parameters) or alternative
ramped harvest control rules can be created with
[make_MP](https://samtool.openmse.com/reference/make_MP.md) and the
available Assess and HCR objects with constant TAC between assessment
years.

## Usage

``` r
SCA_MSY(x, Data, reps = 1, diagnostic = "min")

SCA_75MSY(x, Data, reps = 1, diagnostic = "min")

SCA_4010(x, Data, reps = 1, diagnostic = "min")

DDSS_MSY(x, Data, reps = 1, diagnostic = "min")

DDSS_75MSY(x, Data, reps = 1, diagnostic = "min")

DDSS_4010(x, Data, reps = 1, diagnostic = "min")

SP_MSY(x, Data, reps = 1, diagnostic = "min")

SP_75MSY(x, Data, reps = 1, diagnostic = "min")

SP_4010(x, Data, reps = 1, diagnostic = "min")

SSS_MSY(x, Data, reps = 1, diagnostic = "min")

SSS_75MSY(x, Data, reps = 1, diagnostic = "min")

SSS_4010(x, Data, reps = 1, diagnostic = "min")
```

## Arguments

- x:

  A position in the Data object.

- Data:

  An object of class Data

- reps:

  Numeric, the number of stochastic replicates for the management
  advice.

- diagnostic:

  Character string describing the assessment diagnostic to save, see
  [make_MP](https://samtool.openmse.com/reference/make_MP.md).

## Value

An object of class
[MSEtool::Rec](https://msetool.openmse.com/reference/Rec-class.html)
which contains the management recommendation.

## Functions

- `SCA_MSY()`: A statistical catch-at-age model with a TAC
  recommendation based on fishing at FMSY, and default arguments for
  configuring [SCA](https://samtool.openmse.com/reference/SCA.md).

- `SCA_75MSY()`: An SCA with a TAC recommendation based on fishing at
  75% of FMSY.

- `SCA_4010()`: An SCA with a 40-10 control rule.

- `DDSS_MSY()`: A state-space delay difference model with a TAC
  recommendation based on fishing at FMSY, and default arguments for
  configuring [DD_SS](https://samtool.openmse.com/reference/DD_TMB.md).

- `DDSS_75MSY()`: A state-space delay difference model with a TAC
  recommendation based on fishing at 75% of FMSY.

- `DDSS_4010()`: A state-space delay difference model with a 40-10
  control rule.

- `SP_MSY()`: A surplus production model with a TAC recommendation based
  on fishing at FMSY, and default arguments for configuring
  [SP](https://samtool.openmse.com/reference/SP.md).

- `SP_75MSY()`: A surplus production model with a TAC recommendation
  based on fishing at 75% of FMSY.

- `SP_4010()`: A surplus production model with a 40-10 control rule.

- `SSS_MSY()`: Simple stock synthesis (terminal depletion fixed to 0.4
  in [SSS](https://samtool.openmse.com/reference/SSS.md)) with a TAC
  recommendation based on fishing at FMSY.

- `SSS_75MSY()`: Simple stock synthesis (terminal depletion fixed to
  0.4) with with a TAC recommendation based on fishing at 75% FMSY.

- `SSS_4010()`: Simple stock synthesis (terminal depletion fixed to 0.4)
  with a 40-10 control rule.

## Examples

``` r
MSEtool::avail("MP", package = "SAMtool")
#> → Searching for objects of class "MP" in package "SAMtool"
#>  [1] "DDSS_4010"  "DDSS_75MSY" "DDSS_MSY"   "SCA_4010"   "SCA_75MSY" 
#>  [6] "SCA_MSY"    "SP_4010"    "SP_75MSY"   "SP_MSY"     "SSS_4010"  
#> [11] "SSS_75MSY"  "SSS_MSY"   

# \donttest{
myMSE <- MSEtool::runMSE(MSEtool::testOM, MPs = c("FMSYref", "SCA_4010"))
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
#> → 1/2 Running MSE for: "FMSYref" 
#>   |==                                                | 2 % ~00s            |===                                               | 4 % ~00s            |====                                              | 6 % ~00s            |=====                                             | 8 % ~00s            |======                                            | 10% ~01s            |=======                                           | 12% ~00s            |========                                          | 14% ~00s            |=========                                         | 16% ~00s            |==========                                        | 18% ~00s            |===========                                       | 20% ~00s            |============                                      | 22% ~00s            |=============                                     | 24% ~00s            |==============                                    | 27% ~00s            |===============                                   | 29% ~00s            |================                                  | 31% ~00s            |=================                                 | 33% ~00s            |==================                                | 35% ~00s            |===================                               | 37% ~00s            |====================                              | 39% ~00s            |=====================                             | 41% ~00s            |======================                            | 43% ~00s            |=======================                           | 45% ~00s            |========================                          | 47% ~00s            |=========================                         | 49% ~00s            |==========================                        | 51% ~00s            |===========================                       | 53% ~00s            |============================                      | 55% ~00s            |=============================                     | 57% ~00s            |==============================                    | 59% ~00s            |===============================                   | 61% ~00s            |================================                  | 63% ~00s            |=================================                 | 65% ~00s            |==================================                | 67% ~00s            |===================================               | 69% ~00s            |====================================              | 71% ~00s            |=====================================             | 73% ~00s            |======================================            | 76% ~00s            |=======================================           | 78% ~00s            |========================================          | 80% ~00s            |=========================================         | 82% ~00s            |==========================================        | 84% ~00s            |===========================================       | 86% ~00s            |============================================      | 88% ~00s            |=============================================     | 90% ~00s            |==============================================    | 92% ~00s            |===============================================   | 94% ~00s            |================================================  | 96% ~00s            |================================================= | 98% ~00s            |==================================================| 100% elapsed=01s  
#> 
#> 
#> → 2/2 Running MSE for: "SCA_4010" 
#>   |==                                                | 2 % ~53s            |===                                               | 4 % ~26s            |====                                              | 6 % ~17s            |=====                                             | 8 % ~13s            |======                                            | 10% ~21s            |=======                                           | 12% ~18s            |========                                          | 14% ~15s            |=========                                         | 16% ~13s            |==========                                        | 18% ~17s            |===========                                       | 20% ~15s            |============                                      | 22% ~13s            |=============                                     | 24% ~12s            |==============                                    | 27% ~15s            |===============                                   | 29% ~14s            |================                                  | 31% ~12s            |=================                                 | 33% ~11s            |==================                                | 35% ~14s            |===================                               | 37% ~13s            |====================                              | 39% ~12s            |=====================                             | 41% ~11s            |======================                            | 43% ~13s            |=======================                           | 45% ~12s            |========================                          | 47% ~11s            |=========================                         | 49% ~10s            |==========================                        | 51% ~11s            |===========================                       | 53% ~10s            |============================                      | 55% ~09s            |=============================                     | 57% ~09s            |==============================                    | 59% ~10s            |===============================                   | 61% ~09s            |================================                  | 63% ~08s            |=================================                 | 65% ~07s            |==================================                | 67% ~08s            |===================================               | 69% ~07s            |====================================              | 71% ~07s            |=====================================             | 73% ~06s            |======================================            | 76% ~06s            |=======================================           | 78% ~06s            |========================================          | 80% ~05s            |=========================================         | 82% ~04s            |==========================================        | 84% ~04s            |===========================================       | 86% ~04s            |============================================      | 88% ~03s            |=============================================     | 90% ~03s            |==============================================    | 92% ~02s            |===============================================   | 94% ~02s            |================================================  | 96% ~01s            |================================================= | 98% ~01s            |==================================================| 100% elapsed=30s  
#> 
# }
```
