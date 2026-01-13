# Characterize posterior predictive data

Characterize posterior predictive data

## Usage

``` r
getinds(
  PPD,
  styr,
  res = 6,
  tsd = c("Cat", "Cat", "Cat", "Ind", "ML"),
  stat = c("slp", "AAV", "mu", "slp", "slp")
)
```

## Arguments

- PPD:

  An object of class Data stored in the Misc slot of an MSE object
  following a call of `runMSE(PPD = TRUE)`.

- styr:

  Positive integer, the starting year for calculation of quantities

- res:

  Positive integer, the temporal resolution (chunks - normally years)
  over which to calculate quantities

- tsd:

  Character vector of names of types of data: Cat = catch, Ind =
  relative abundance index, ML = mean length in catches

- stat:

  Character vector of types of quantity to be calculated: slp =
  slope(log(x)), AAV = average annual variability, mu = mean(log(x))

## Value

A 3D array of results (type of data/stat (e.g. mean catches),time period
(chunk), simulation)

## References

Carruthers and Hordyk 2018

## Author

T. Carruthers
