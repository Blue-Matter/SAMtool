# retrospective_AM (retrospective of Assessment model in MSE)

Plots the true retrospective of an assessment model during the
closed-loop simulation. A series of time series estimates of SSB, F, and
VB are plotted over the course of the MSE are plotted against the
operating model (true) values (in black).

## Usage

``` r
retrospective_AM(MSE, MP, sim = 1, plot_legend = FALSE)
```

## Arguments

- MSE:

  An object of class
  [MSEtool::MSE](https://msetool.openmse.com/reference/MSE-class.html).

- MP:

  Character. The name of the management procedure created by
  [`make_MP()`](https://samtool.openmse.com/reference/make_MP.md)
  containing the assessment model.

- sim:

  Integer between 1 and MSE@nsim. The simulation number for which the
  retrospectives will be plotted.

- plot_legend:

  Logical. Whether to plot legend to reference year of assessment in the
  MSE.

## Value

A series of figures for SSB, depletion, fishing mortality, and
vulnerable biomass (VB) estimated in the MP over the course of the
closed-loop simulation against the values generated in the operating
model (both historical and projected).

## Details

For assessment models that utilize annual exploitation rates (u), the
instantaneous fishing mortality rates are obtained as F = -log(1 - u).

## Note

This function only plots retrospectives from a single simulation in the
MSE. Results from one figure may not be indicative of general assessment
behavior and performance overall.

## See also

[diagnostic](https://samtool.openmse.com/reference/diagnostic.md)

## Author

Q. Huynh

## Examples

``` r
# \donttest{
SP_40_10 <- make_MP(SP, HCR_MSY, diagnostic = "full")
OM <- MSEtool::testOM; OM@proyears <- 20
myMSE <- MSEtool::runMSE(OM = OM, MPs = "SP_40_10")
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
#> → 1/1 Running MSE for: "SP_40_10" 
#>   |===                                               | 5 % ~02s            |======                                            | 11% ~01s            |========                                          | 16% ~01s            |===========                                       | 21% ~01s            |==============                                    | 26% ~01s            |================                                  | 32% ~01s            |===================                               | 37% ~01s            |======================                            | 42% ~01s            |========================                          | 47% ~01s            |===========================                       | 53% ~00s            |=============================                     | 58% ~00s            |================================                  | 63% ~00s            |===================================               | 68% ~00s            |=====================================             | 74% ~00s            |========================================          | 79% ~00s            |===========================================       | 84% ~00s            |=============================================     | 89% ~00s            |================================================  | 95% ~00s            |==================================================| 100% elapsed=01s  
#> 
retrospective_AM(myMSE, MP = "SP_40_10", sim = 1)


# How to get all the estimates
library(dplyr)
assess_estimates <- lapply(1:myMSE@nMPs, function(m) {
  lapply(1:myMSE@nsim, function(x) {
    report <- myMSE@PPD[[m]]@Misc[[x]]$Assessment_report
    if (is.null(report)) {
      return(data.frame())
    } else {
      mutate(report, MP = myMSE@MPs[m], Simulation = x)
    }
  }) %>% bind_rows()
}) %>% bind_rows()
# }
```
