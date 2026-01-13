# Projections for assessment models

This function takes an assessment model and runs a stochastic projection
based on future F or catch.

## Usage

``` r
projection(
  Assessment,
  constrain = c("F", "Catch"),
  Ftarget,
  Catch,
  p_years = 50,
  p_sim = 200,
  obs_error,
  process_error,
  max_F = 3,
  seed = 499
)
```

## Arguments

- Assessment:

  An object of class
  [Assessment](https://samtool.openmse.com/reference/Assessment-class.md).

- constrain:

  Whether to project on future F or catch. By default, projects on F.

- Ftarget:

  The projection F, either of length 1 for constant F for the entirety
  of the projection or length p_years.

- Catch:

  The projection catch, either of length 1 for constant catch for the
  entirety of the projection or length p_years.

- p_years:

  Integer for the number of projection years.

- p_sim:

  Integer for the number of simulations for the projection.

- obs_error:

  A list of length two. In the first entry, a vector of length nsurvey
  giving the standard deviations of each future index, or alternatively
  an array of dimension p_sim, p_years, and nsurvey giving the deviates.
  The second entry is the standard deviation of the projected catch.
  Alternatively, a matrix of simulation and year-specific error
  structure for the catch (p_sim rows and p_year columns; a matrix of
  ones indicates perfect data).

- process_error:

  Numeric, standard deviation for process error (e.g., recruitment or
  biomass deviates). If `NULL`, uses values from assessment model.
  Alternatively, a matrix of simulation and year-specific recruitment
  deviates (p_sim rows and p_year columns, a matrix of ones indicates no
  recruitment deviates).

- max_F:

  The maximum allowable F if the projection is constrained on catch.

- seed:

  An integer to set the seed for the sampling observation and process
  error deviates.

## Value

An object of class
[project](https://samtool.openmse.com/reference/project-class.md) that
contains future predicted values of F, catch, biomass, recruitment, etc.

## Examples

``` r
# \donttest{
myAssess <- SP(Data = swordfish)
do_projection <- projection(myAssess, Ftarget = myAssess@FMSY)
# }
```
