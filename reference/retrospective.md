# Retrospective analysis of assessment models

Perform a retrospective analysis, successive removals of most recent
years of data to evaluate resulting parameter estimates.

## Usage

``` r
retrospective(x, ...)

# S4 method for class 'Assessment'
retrospective(x, nyr = 5, figure = TRUE)

# S4 method for class 'RCModel'
retrospective(x, nyr = 5, figure = TRUE)
```

## Arguments

- x:

  An S4 object of class
  [Assessment](https://samtool.openmse.com/reference/Assessment-class.md)
  or [RCModel](https://samtool.openmse.com/reference/RCModel-class.md).

- ...:

  More arguments.

- nyr:

  The maximum number of years to remove for the retrospective analysis.

- figure:

  Indicates whether plots will be drawn.

## Value

A list with an array of model output and of model estimates from the
retrospective analysis.

Figures showing the time series of biomass and exploitation and
parameter estimates with successive number of years removed. For a
variety of time series output (SSB, recruitment, etc.) and estimates
(R0, steepness, etc.), also returns a matrix of Mohn's rho (Mohn 1999).

## References

Mohn, R. 1999. The retrospective problem in sequential population
analysis: an investigation using cod fishery and simulated data. ICES
Journal of Marine Science 56:473-488.

## Author

Q. Huynh

## Examples

``` r
# \donttest{
output <- SP(Data = swordfish)
get_retro <- retrospective(output, nyr = 5, figure = FALSE)
# }
```
