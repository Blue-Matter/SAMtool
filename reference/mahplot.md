# Plot statistical power of the indicator with increasing time blocks

Plot statistical power of the indicator with increasing time blocks

## Usage

``` r
mahplot(outlist, res = 6, maxups = 5, MPs)
```

## Arguments

- outlist:

  A list object produced by the function
  [PRBcalc](https://samtool.openmse.com/reference/PRBcalc.md)

- res:

  Integer, the resolution (time blocking) for the calculation of PPD

- maxups:

  Integer, the maximum number of update time blocks to plot

- MPs:

  Character vector of MP names

## Value

Density plots of Mahalanobis distance.

## References

Carruthers and Hordyk 2018

## Author

T. Carruthers
