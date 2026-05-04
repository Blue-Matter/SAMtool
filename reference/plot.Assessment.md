# Plot Assessment object

Produces HTML file (via markdown) figures of parameter estimates and
output from an
[Assessment](https://samtool.openmse.com/reference/Assessment-class.md)
object.

## Usage

``` r
# S4 method for class 'Assessment,missing'
plot(
  x,
  filename = paste0("report_", x@Model),
  dir = tempdir(),
  ret_yr = 0L,
  open_file = TRUE,
  quiet = TRUE,
  render_args = list(),
  ...
)

# S4 method for class 'Assessment,retro'
plot(
  x,
  y,
  filename = paste0("report_", x@Model),
  dir = tempdir(),
  open_file = TRUE,
  quiet = TRUE,
  render_args = list(),
  ...
)
```

## Arguments

- x:

  An object of class
  [Assessment](https://samtool.openmse.com/reference/Assessment-class.md).

- filename:

  Character string for the name of the markdown and HTML files.

- dir:

  The directory in which the markdown and HTML files will be saved.

- ret_yr:

  If greater than zero, then a retrospective analysis will be performed
  and results will be reported. The integer here corresponds to the
  number of peels (the maximum number of terminal years for which the
  data are removed).

- open_file:

  Logical, whether the HTML document is opened after it is rendered.

- quiet:

  Logical, whether to silence the markdown rendering function.

- render_args:

  Arguments to pass to
  [render](https://pkgs.rstudio.com/rmarkdown/reference/render.html).

- ...:

  Other arguments.

- y:

  An object of class
  [retro](https://samtool.openmse.com/reference/retro-class.md).

## Value

Returns invisibly the output from
[render](https://pkgs.rstudio.com/rmarkdown/reference/render.html).

## See also

[retrospective](https://samtool.openmse.com/reference/retrospective.md)

## Examples

``` r
output <- DD_TMB(Data = Simulation_1)

# \donttest{
plot(output)
#> ℹ Writing markdown file: /tmp/RtmpuYmaxh/report_DD_TMB.Rmd
#> ℹ Rendering markdown file: /tmp/RtmpuYmaxh/report_DD_TMB.Rmd
#> ℹ See help(plot.Assessment) to adjust report and file directory.
#> ✔ Rendered file: /tmp/RtmpuYmaxh/report_DD_TMB.html
# }
```
