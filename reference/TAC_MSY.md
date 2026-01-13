# Calculate MSY-based TAC from Assessment object

A function to calculate the total allowable catch (TAC). Based on the
MSY (maximum sustainable yield) principle, the TAC is the product of
either UMSY or FMSY and the available biomass, i.e. vulnerable biomass,
in terminal year.

## Usage

``` r
TAC_MSY(Assessment, reps, MSY_frac = 1)
```

## Arguments

- Assessment:

  An Assessment object with estimates of UMSY or FMSY and terminal year
  vulnerable biomass.

- reps:

  The number of stochastic draws of UMSY or FMSY.

- MSY_frac:

  The fraction of FMSY or UMSY for calculating the TAC (e.g. MSY_frac =
  0.75 fishes at 75% of FMSY).

## Value

A vector of length `reps` of stochastic samples of TAC recommendation.
Returns NA's if missing either UMSY/FMSY or vulnerable biomass.

## Note

`calculate_TAC` is deprecated as of version 1.2 in favor of `TAC_MSY`
because the latter has a more informative name.

## See also

[HCR_MSY](https://samtool.openmse.com/reference/HCR_MSY.md)
[HCR40_10](https://samtool.openmse.com/reference/HCR_ramp.md)
[HCR60_20](https://samtool.openmse.com/reference/HCR_ramp.md)
