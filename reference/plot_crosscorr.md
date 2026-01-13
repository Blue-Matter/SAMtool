# Produce a cross-correlation plot of the derived data arising from getinds(MSE_object)

Produce a cross-correlation plot of the derived data arising from
getinds(MSE_object)

## Usage

``` r
plot_crosscorr(
  indPPD,
  indData,
  pp = 1,
  dnam = c("CS", "CV", "CM", "IS", "MLS"),
  res = 1
)
```

## Arguments

- indPPD:

  A 3D array of results arising from running getind on an MSE of the
  Null operating model (type of data/stat (e.g. mean catches),time
  period (chunk), simulation)

- indData:

  A 3D array of results arising from running getind on an MSE of the
  Alternative operating model (type of data/stat (e.g. mean
  catches),time period (chunk), simulation)

- pp:

  Positive integer, the number of time chunks (blocks of years normally,
  second dimension of indPPD and indData) to produce the plot for.

- dnam:

  A character vector of names of the data for plotting purposes (as long
  as dimension 1 of indPPD and indData).

- res:

  The size of the temporal blocking that created indPPD and indData -
  this is just used for labelling purposes

## Value

A cross-correlation plot (ndata-1) x (ndata-1)

## References

Carruthers and Hordyk 2018

## Author

T. Carruthers
