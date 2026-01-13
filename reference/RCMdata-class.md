# Class-`RCMdata`

An S4 class for the data inputs into
[RCM](https://samtool.openmse.com/reference/RCM.md).

## Slots

- `Chist`:

  Either a vector of historical catch, should be of length `OM@nyears`,
  or if there are multiple fleets, a matrix of `OM@nyears` rows and
  `nfleet` columns. Ideally, the first year of the catch series
  represents unfished conditions (see also slot `C_eq`).

- `C_sd`:

  Same dimension as `Chist`. Lognormal distribution standard deviations
  (by year and fleet) for the catches in `Chist`. If not provided, the
  default is 0.01. Not used if `RCM(condition = "catch2")`.

- `Ehist`:

  A vector of historical effort, should be of length `OM@nyears`, or if
  there are multiple fleets: a matrix of `OM@nyears` rows and `nfleet`
  columns. See also slot `E_eq`).

- `C_wt`:

  Optional weight at age for the catch `Chist`. Array with dimension
  `[OM@nyears+1, OM@maxage+1, nfleet]`.

- `CAA`:

  Fishery age composition matrix with `nyears` rows and `OM@maxage+1`
  columns, or if multiple fleets: an array with dimension:
  `nyears, OM@maxage+1, nfleet`. Enter `NA` for years without any data.
  Raw numbers will be converted to annual proportions (see slot
  `CAA_ESS` for sample sizes).

- `CAA_ESS`:

  Annual sample size (for the multinomial distribution) of the fishery
  age comps. A vector of length `OM@nyears`, or if there are multiple
  fleets: a matrix of `OM@nyears` rows and `nfleet` columns. Enter zero
  for years without observations. An annual cap to the ESS, e.g., 50,
  can be calculated with something like:
  `pmin(apply(CAA, c(1, 3), sum, na.rm = TRUE), 50)`. By default,

- `CAL`:

  Fishery length composition matrix with `nyears` rows and `n_bin`
  columns (indexing the length bin), or if multiple fleets: an array
  with dimension: `nyears, n_bin, nfleets`. Enter `NA` for years without
  any data. Raw numbers will be converted to annual proportions (see
  slot `CAL_ESS` for sample sizes).

- `CAL_ESS`:

  Annual sample size (for the multinomial distribution) of the fishery
  length comps. Same dimension as `CAA_ESS`.

- `length_bin`:

  - A vector (length `n_bin`) for the midpoints of the length bins for
    `CAL` and `IAL`, as well as the population model, if all bin widths
    are equal in size. If length bins are unequal in width, then provide
    a vector of the boundaries of the length bins (vector of length
    `n_bin + 1`).

- `MS`:

  Mean mean size (either mean length or mean weight) observations from
  the fishery. Same dimension as `Chist`. Generally, mean lengths should
  not be used alongside `CAL`, unless mean length and length comps are
  independently sampled.

- `MS_type`:

  A character (either `"length"` (default) or `"weight"`) to denote the
  type of mean size data.

- `MS_cv`:

  The coefficient of variation of the observed mean size. If there are
  multiple fleets, a vector of length `nfleet`. Default is 0.2.

- `Index`:

  Index of abundance. Enter `NA` for missing values. A vector length
  `OM@nyears`, or if there are multiple surveys: a matrix of `OM@nyears`
  rows and `nsurvey` columns.

- `I_sd`:

  A vector or matrix of standard deviations (lognormal distribution) for
  the indices corresponding to the entries in `Index`. Same dimension as
  `Index`. If not provided, this function will use values from
  `OM@Iobs`.

- `I_wt`:

  Optional weight at age for the index `Index`. Array with dimension
  `[OM@nyears, OM@maxage+1, nsurvey]`.

- `IAA`:

  Index age composition data, an array of dimension
  `nyears, maxage+1, nsurvey`. Raw numbers will be converted to annual
  proportions (see `IAA_ESS` for sample sizes).

- `IAA_ESS`:

  Annual sample size (for the multinomial distribution) of the index age
  comps. A vector of length `OM@nyears`. If there are multiple indices:
  a matrix of `OM@nyears` rows and `nsurvey` columns.

- `IAL`:

  Index length composition data, an array of dimension
  `nyears, n_bin, nsurvey`. Raw numbers will be converted to annual
  proportions (see slot `IAL_ESS` to enter sample sizes).

- `IAL_ESS`:

  Annual sample size (for the multinomial distribution) of the index
  length comps. Same dimension as `IAA_ESS`.

- `C_eq`:

  Vector of length `nfleet` for the equilibrium catch for each fleet in
  `Chist` prior to the first year of the operating model. Zero (default)
  implies unfished conditions in year one. Otherwise, this is used to
  estimate depletion in the first year of the data. Alternatively, if
  one has a full CAA matrix, one could instead estimate "artificial" rec
  devs to generate the initial numbers-at-age (and hence initial
  depletion) in the first year of the model (see additional arguments in
  [RCM](https://samtool.openmse.com/reference/RCM.md)).

- `C_eq_sd`:

  - A vector of standard deviations (lognormal distribution) for the
    equilibrium catches in `C_eq`. Same dimension as `C_eq`. If not
    provided, the default is 0.01. Only used if
    `RCM(condition = "catch")`.

- `E_eq`:

  The equilibrium effort for each fleet in `Ehist` prior to the first
  year of the operating model. Zero (default) implies unfished
  conditions in year one. Otherwise, this is used to estimate depletion
  in the first year of the data.

- `abs_I`:

  An integer vector length `nsurvey` to indicate which indices are in
  absolute magnitude. Use `1` to set `q = 1`, otherwise use 0 (default)
  to estimate q.

- `I_units`:

  An integer vector of length `nsurvey` to indicate whether indices are
  biomass based (1) or abundance-based (0). By default, all are
  biomass-based.

- `I_delta`:

  A vector of length `nsurvey` to indicate the timing of the indices
  within each time step (0-1, for example 0.5 is the midpoint of the
  year). By default, zero is used. Can also be a matrix by
  `nyears, nsurvey`. Use -1 if the survey operates continuously, the
  availability would be `N * (1 - exp(-Z))/Z`.

- `age_error`:

  A square matrix of `maxage + 1` rows and columns to specify ageing
  error. The `aa`-th column assigns a proportion of animals of true age
  `aa` to observed age `a` in the `a`-th row. Thus, all rows should sum
  to 1. Default is an identity matrix (no ageing error).

- `sel_block`:

  For time-varying fleet selectivity (in time blocks), a integer matrix
  of `nyears` rows and `nfleet` columns to assign a selectivity function
  to a fleet for certain years. By default, constant selectivity for
  each individual fleet. See the
  [selectivity](https://openmse.com/tutorial-rcm-select/) article for
  more details.

- `Misc`:

  A list of miscellaneous inputs. Used internally.

## See also

[RCM](https://samtool.openmse.com/reference/RCM.md)

## Author

Q. Huynh
