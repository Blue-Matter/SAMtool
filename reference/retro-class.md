# Class-`retro`

An S4 class that contains output from
[retrospective](https://samtool.openmse.com/reference/retrospective.md).

## Slots

- `Model`:

  Name of the assessment model.

- `Name`:

  Name of Data object.

- `TS_var`:

  Character vector of time series variables, e.g. recruitment, biomass,
  from the assessment.

- `TS`:

  An array of time series assessment output of dimension, indexed by:
  peel (the number of terminal years removed from the base assessment),
  years, and variables (corresponding to `TS_var`).

- `Est_var`:

  Character vector of estimated parameters, e.g. R0, steepness, in the
  assessment.

- `Est`:

  An array for estimated parameters of dimension, indexed by: peel,
  variables (corresponding to `Est_var`), and value (length 2 for
  estimate and standard error).

## See also

[plot.retro](https://samtool.openmse.com/reference/plot.retro.md)
[summary.retro](https://samtool.openmse.com/reference/plot.retro.md)
[plot.Assessment](https://samtool.openmse.com/reference/plot.Assessment.md)

## Author

Q. Huynh
