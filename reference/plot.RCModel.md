# Plot RCM scope output

Produces HTML file (via markdown) figures of parameter estimates and
output from an
[Assessment](https://samtool.openmse.com/reference/Assessment-class.md)
object. Plots histograms of operating model parameters that are updated
by the RCM scoping function, as well as diagnostic plots for the fits to
the RCM for each simulation. `compare_RCM` plots a short report that
compares output from multiple RCM objects, assuming the same model
structure, i.e., identical matrix and array dimensions among models, but
different data weightings, data omissions, etc.

## Usage

``` r
# S4 method for class 'RCModel,missing'
plot(
  x,
  compare = FALSE,
  filename = "RCM",
  dir = tempdir(),
  sims = 1:x@OM@nsim,
  Year = NULL,
  Age = NULL,
  f_name = NULL,
  s_name = NULL,
  MSY_ref = c(0.5, 1),
  bubble_adj = 1.5,
  scenario = list(),
  title = NULL,
  open_file = TRUE,
  quiet = TRUE,
  render_args,
  ...
)

compare_RCM(
  ...,
  compare = FALSE,
  filename = "compare_RCM",
  dir = tempdir(),
  Year = NULL,
  Age = NULL,
  f_name = NULL,
  s_name = NULL,
  MSY_ref = c(0.5, 1),
  bubble_adj = 1.5,
  scenario = list(),
  title = NULL,
  open_file = TRUE,
  quiet = TRUE,
  render_args
)
```

## Arguments

- x:

  An object of class
  [RCModel](https://samtool.openmse.com/reference/RCModel-class.md)
  (output from [RCM](https://samtool.openmse.com/reference/RCM.md)).

- compare:

  Logical, if TRUE, the function will run `runMSE` to compare the
  historical period of the operating model and the RCM output.

- filename:

  Character string for the name of the markdown and HTML files.

- dir:

  The directory in which the markdown and HTML files will be saved.

- sims:

  A logical vector of length `x@OM@nsim` or a numeric vector indicating
  which simulations to keep.

- Year:

  Optional, a vector of years for the historical period for plotting.
  Useful if seasonal time steps are used.

- Age:

  Optional, a vector of ages for plotting. Useful if seasonal time steps
  are used.

- f_name:

  Character vector for fleet names.

- s_name:

  Character vector for survey names.

- MSY_ref:

  A numeric vector for reference horizontal lines for B/BMSY plots.

- bubble_adj:

  A number to adjust the size of bubble plots (for residuals of age and
  length comps).

- scenario:

  Optional, a named list to label each simulation in the RCM for
  plotting, e.g.:
  `list(names = c("low M", "high M"), col = c("blue", "red"))`.

- title:

  Optional character string for an alternative title for the markdown
  report.

- open_file:

  Logical, whether the HTML document is opened after it is rendered.

- quiet:

  Logical, whether to silence the markdown rendering function.

- render_args:

  A list of other arguments to pass to
  [render](https://pkgs.rstudio.com/rmarkdown/reference/render.html).

- ...:

  For `compare_RCM`, multiple RCM objects for comparison.

## Value

Returns invisibly the output from
[render](https://pkgs.rstudio.com/rmarkdown/reference/render.html).

## See also

[RCModel](https://samtool.openmse.com/reference/RCModel-class.md)
[RCM](https://samtool.openmse.com/reference/RCM.md)
