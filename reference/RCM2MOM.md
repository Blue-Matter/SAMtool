# Convert RCM to a multi-fleet operating model (MOM)

The RCM (Rapid Conditioning Model) returns a single-fleet operating
model, implying constant effort among fleets for projections. Here, we
convert the single-fleet OM to a multi-fleet OM, preserving the multiple
fleet structure used in the conditioning model for projections. This
allows for testing management procedures that explicitly specify fleet
allocation in the management advice.

## Usage

``` r
RCM2MOM(RCModel)
```

## Arguments

- RCModel:

  Output from [RCM](https://samtool.openmse.com/reference/RCM.md), a
  class
  [RCModel](https://samtool.openmse.com/reference/RCModel-class.md)
  object.

## Value

A class
[MSEtool::MOM](https://msetool.openmse.com/reference/MOM-class.html)
object.

## Author

Q. Huynh

## Examples

``` r
# \donttest{
data(pcod) 
mat_ogive <- pcod$OM@cpars$Mat_age[1, , 1]
OM <- MSEtool::SubCpars(pcod$OM, 1:3)
#> ✔ Removing simulations:  4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26 27 28 29 30 31 32 33 34 35 36 37 38 39 40 41 42 43 44 45 46 47 48
#> ✔ Set OM@nsim =  3
out <- RCM(OM = pcod$OM, data = pcod$data, 
           condition = "catch", mean_fit = TRUE,
           selectivity = "free", s_selectivity = rep("SSB", ncol(pcod$data@Index)),
           start = list(vul_par = matrix(mat_ogive, length(mat_ogive), 1)),
           map = list(vul_par = matrix(NA, length(mat_ogive), 1),
                      log_early_rec_dev = rep(1, pcod$OM@maxage)),
           prior = pcod$prior)
#> ℹ Checking data...
#> ✔ 1 fleet(s) detected.
#> ✔ RCM is conditioned on:
#> ✔ Fleet 1: catch
#> ✔ 65 years of data detected.
#> ✔ First year in model: 1956
#> ✔ Last year in model: 2020
#> ✔ 5 survey(s) detected.
#> ✔ Checking OM and getting biological parameters...
#> ✔ Mean weight data found.
#> ✔ Passing user arguments (LWT, map, start, prior, etc.) to RCMdata@Misc..
#> ✔ Maximum F in RCM will be 3. OM@maxF is also updated.
#>   
#> ℹ No fishery length or age compositions were provided. Selectivity is fixed to values from OM.
#>   
#> ℹ Fishery selectivity setup:
#> ℹ Fleet 1: individual parameters at age (free)
#>   
#> ℹ Index selectivity setup:
#> ℹ Index 1: spawning biomass
#> ℹ Index 2: spawning biomass
#> ℹ Index 3: spawning biomass
#> ℹ Index 4: spawning biomass
#> ℹ Index 5: spawning biomass
#>   
#> ✔ Beverton-Holt stock-recruitment relationship used.
#> ✔ Prior for q found.
#> ℹ Fitting model (48 simulations) ...
#> ℹ Generating additional model fit from mean values of parameters in the operating model...
#> ✔ Updating operating model:
#>   
#> ✔ Range of unfished age-0 recruitment (OM@cpars$R0): 6383.46 - 13782.9
#> ✔ Range of initial spawning depletion: 0.44 - 1.45
#> ✔ Range of spawning depletion (OM@cpars$D): 0.14 - 0.45
#> ✔ Historical F set with OM@cpars$Find and OM@cpars$qs.
#> ✔ Annual selectivity at age set in OM@cpars$V. Projection period uses selectivity of last historical year.
#> ✔ RCMdata length bins will be added to OM.
#> ✔ Recruitment standard deviation set in OM@cpars$Perr: 0.8 - 0.8
#> ✔ Historical recruitment deviations set in OM@cpars$Perr_y.
#> ✔ Range of recruitment autocorrelation OM@AC: 0.21 - 0.31
#> ✔ Future recruitment deviations in OM@cpars$Perr_y sampled with autocorrelation.
#> ✔ Growth, maturity, natural mortality, and stock recruit parameters from RCM are set in OM@cpars.
#>   
#> ℹ Adding some RCMdata inputs into OM@cpars$Data:
#>   
#> ✔ Historical catch data added to OM@cpars$Data@Cat.
#> ✔ Historical indices added to OM@cpars$Data@AddInd.
#> ✔ Complete.
MOM <- RCM2MOM(out)
# }
```
