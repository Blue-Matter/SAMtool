# Class-`project`

An S4 class for the output from
[projection](https://samtool.openmse.com/reference/projection.md).

## Slots

- `Model`:

  Name of the assessment model.

- `Name`:

  Name of Data object.

- `FMort`:

  A matrix of fishing mortality over `p_sim` rows and `p_years` columns.

- `B`:

  An matrix of biomass with `p_sim` rows and `p_years` columns.

- `SSB`:

  A matrix of spawning biomass with `p_sim` rows and `p_years` columns.

- `VB`:

  A matrix of vulnerable biomass with `p_sim` rows and `p_years`
  columns.

- `R`:

  A matrix of recruitment over `p_sim` rows and `p_years` columns.

- `N`:

  A matrix of abundance over `p_sim` rows and `p_years` columns.

- `Catch`:

  A matrix of simulated observed catch over `p_sim` rows and `p_years`
  columns.

- `Index`:

  An array of simulated observed index of dimension
  `c(p_sim, p_years, nsurvey)`.

- `C_at_age`:

  An array for catch-at-age with dimension `c(p_sim, p_years, n_age)`.

## See also

[projection](https://samtool.openmse.com/reference/projection.md)

## Author

Q. Huynh
