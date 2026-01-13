# Class-`prof`

An S4 class that contains output from
[profile](https://samtool.openmse.com/reference/profile.md).

## Slots

- `Model`:

  Name of the assessment model.

- `Name`:

  Name of Data object.

- `Par`:

  Character vector of parameters that were profiled.

- `MLE`:

  Numeric vector of the estimated values of the parameters
  (corresponding to `Par`) from the assessment.

- `grid`:

  A data.frame of the change in negative log-likelihood (`nll`) based on
  the profile of the parameters.

## See also

[plot.prof](https://samtool.openmse.com/reference/plot.prof.md)
[profile](https://samtool.openmse.com/reference/profile.md)

## Author

Q. Huynh
