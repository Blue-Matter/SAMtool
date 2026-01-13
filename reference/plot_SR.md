# Plot stock-recruitment function

Plot stock-recruitment (with recruitment deviations if estimated).

## Usage

``` r
plot_SR(
  Spawners,
  expectedR,
  R0 = NULL,
  S0 = NULL,
  rec_dev = NULL,
  trajectory = FALSE,
  y_zoom = NULL,
  ylab = "Recruitment"
)
```

## Arguments

- Spawners:

  A vector of the number of the spawners (x-axis).

- expectedR:

  A vector of the expected recruitment (from the stock-recruit function)
  corresponding to values of `Spawners`.

- R0:

  Virgin recruitment.

- S0:

  Virgin spawners.

- rec_dev:

  If recruitment deviations are estimated, a vector of estimated
  recruitment (in normal space) corresponding to values of `Spawners`.

- trajectory:

  Indicates whether arrows will be drawn showing the trajectory of
  spawners and recruitment deviations over time.

- y_zoom:

  If recruitment deviations are plotted, the y-axis limit relative to
  maximum expected recruitment `expectedR`. If `NULL`, all recruitment
  values are plotted.

- ylab:

  Character string for label on y-axis.

## Value

A stock-recruit plot

## Author

Q. Huynh
