# Plot residuals

Plots figure of residuals (or any time series with predicted mean of
zero).

## Usage

``` r
plot_residuals(
  Year,
  res,
  res_sd = NULL,
  res_sd_CI = 0.95,
  res_upper = NULL,
  res_lower = NULL,
  res_ind_blue = NULL,
  draw_zero = TRUE,
  zero_linetype = 2,
  label = "Residual"
)
```

## Arguments

- Year:

  A vector of years for the data.

- res:

  A vector of residuals.

- res_sd:

  A vector of year specific standard deviation for `res`.

- res_sd_CI:

  The confidence interval for the error bars based for `res_sd`.

- res_upper:

  A vector of year-specific upper bounds for the error bars of the
  residual (in lieu of argument `res_CV`).

- res_lower:

  A vector of year-specific lower bounds for the error bars of the
  residual (in lieu of argument `res_CV`).

- res_ind_blue:

  Indices of `obs` for which the plotted residuals and error bars will
  be blue.

- draw_zero:

  Indicates whether a horizontal line should be drawn at zero.

- zero_linetype:

  Passes argument `lty` (e.g. solid line = 1, dotted = 2) to
  `draw_zero`.

- label:

  Character string that describes the data to label the y-axis.

## Value

A plot of model residuals by year (optionally, with error bars).

## See also

[`plot_timeseries()`](https://samtool.openmse.com/reference/plot_timeseries.md)

## Author

Q. Huynh
