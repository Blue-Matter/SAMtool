# Summary of Assessment object

Returns a summary of parameter estimates and output from an
[Assessment](https://samtool.openmse.com/reference/Assessment-class.md)
object.

## Usage

``` r
# S4 method for class 'Assessment'
summary(object)
```

## Arguments

- object:

  An object of class
  [Assessment](https://samtool.openmse.com/reference/Assessment-class.md)

## Value

A list of parameters.

## Examples

``` r
output <- DD_TMB(Data = MSEtool::SimulatedData)
summary(output)
#> $model
#> [1] "Delay-Difference"
#> 
#> $current_status
#>            Value
#> F/FMSY 0.5416090
#> B/BMSY 1.5522568
#> B/B0   0.3921528
#> 
#> $input_parameters
#>            Value                         Description
#> alpha  9.3839612              alpha = Winf * (1-rho)
#> rho    0.8460155 rho = (W_k+2 - Winf)/(W_k+1 - Winf)
#> k      4.0000000       Age of knife-edge selectivity
#> w_k   21.3874412                     Weight at age k
#> h      0.7341829             Stock-recruit steepness
#> M      0.3835144                   Natural mortality
#> 
#> $derived_quantities
#>                Value                     Description
#> B0      1.520520e+04                Unfished biomass
#> N0      4.875238e+02              Unfished abundance
#> MSY     2.059243e+03 Maximum sustainable yield (MSY)
#> FMSY    9.834573e-01        Fishing mortality at MSY
#> BMSY    3.841350e+03                  Biomass at MSY
#> BMSY/B0 2.526340e-01                Depletion at MSY
#> 
#> $model_estimates
#>                Estimate   Std. Error     Gradient          CV
#> R0         1.552947e+02 1.855504e+00           NA 0.011948276
#> h          7.341829e-01 0.000000e+00           NA 0.000000000
#> q          1.096571e-04 3.111379e-06           NA 0.028373713
#> sigma      3.890664e-01 5.502230e-03           NA 0.014142136
#> R0x       -2.391703e+00 1.194828e-02 8.796344e-05 0.004995719
#> log_sigma -9.440052e-01 1.414214e-02 7.732051e-05 0.014980994
#> 
#> $log_likelihood
#>           Neg.LL
#> Total   1187.333
#> Index_1 1187.333
#> Dev        0.000
#> Prior      0.000
#> Penalty    0.000
#> 
```
