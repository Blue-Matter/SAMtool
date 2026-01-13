# Segmented harvest control rules

A linear segmented output control rule where the target F (used for the
TAC recommendation) is a function of an operational control point (OCP)
such as spawning depletion or spawning biomass. The segments of the HCR
are specified by arguments `OCP` and `relF`. Beyond the range of `OCP`,
the response will be flat.
[HCR_ramp](https://samtool.openmse.com/reference/HCR_ramp.md) uses
`HCR_segment` with two control points.

## Usage

``` r
HCR_segment(
  Assessment,
  reps = 1,
  OCP_type = c("SSB_SSB0", "SSB_SSBMSY", "SSB_dSSB0", "F_FMSY", "F_F01", "F_FSPR"),
  Ftarget_type = c("FMSY", "F01", "Fmax", "FSPR", "abs"),
  OCP = c(0.1, 0.4),
  relF = c(0, 1),
  SPR_OCP,
  SPR_targ,
  ...
)
```

## Arguments

- Assessment:

  An object of class
  [Assessment](https://samtool.openmse.com/reference/Assessment-class.md)
  with estimates of FMSY or UMSY, vulnerable biomass, and spawning
  biomass depletion in terminal year.

- reps:

  The number of stochastic samples of the TAC recommendation.

- OCP_type:

  The type of operational control points (OCPs) for the harvest control
  rule used to determine the reduction in F. See below.

- Ftarget_type:

  The type of F used for the target fishing mortality rate. See below.

- OCP:

  Numeric vector of operational control points for the HCR (in
  increasing order).

- relF:

  Numeric vector of Ftarget corresponding to the values in `OCP`.

- SPR_OCP:

  The value of spawning potential ratio for the OCP if
  `OCP_type = "F_FSPR"`. By default, 0.4 (F40%).

- SPR_targ:

  The target value of spawning potential ratio if
  `Ftarget_type = "FSPR"`. By default, 0.4 (F40%).

- ...:

  Miscellaneous arguments.

## Value

An object of class
[MSEtool::Rec](https://msetool.openmse.com/reference/Rec-class.html)
with the TAC recommendation.

## Details

The catch advice is calculated using the catch equation of the
corresponding assessment. See `Assessment@forecast$catch_eq`, a function
that returns the catch advice for a specified `Ftarget`.

*Operational control points (OCP_type)*

The following are the available options for harvest control rule inputs,
and the source of those values in the
[Assessment](https://samtool.openmse.com/reference/Assessment-class.md)
object:

- **Default** `"SSB_SSB0"`: Spawning depletion. Uses the last value in
  `Assessment@SSB_SSB0` vector.

- `"SSB_SSBMSY"`: Spawning biomass relative to MSY. Uses the last value
  in `Assessment@SSB_SSBMSY` vector.

- `"SSB_dSSB0"`: Dynamic depletion (SSB relative to the historical
  reconstructed biomass with F = 0). Uses the last value in
  `Assessment@SSB/Assessment@TMB_report$dynamic_SSB0`.

- `"F_FMSY"`: Fishing mortality relative to MSY. Uses the last value in
  `Assessment@F_FMSY`.

- `"F_F01"`: Fishing mortality relative to F_0.1 (yield per recruit),
  calculated from the data frame in
  `Assessment@forecast[["per_recruit"]]`.

- `"F_FSPR"`: Fishing mortality relative to F_SPR% (the F that produces
  the spawning potential ratio specified in `"SPR_OCP"`, calculated from
  the data frame in `Assessment@forecast[["per_recruit"]]`.

*Fishing mortality target (Ftarget_type)*

The type of F for which the corresponding catch is calculated in the HCR
is specified here. The source of those values in the
[Assessment](https://samtool.openmse.com/reference/Assessment-class.md)
object is specified:

- **Default** `"FMSY"`: Fishing mortality relative to MSY. Uses the
  value in `Assessment@FMSY`.

- `"F01"`: Fishing mortality relative to F_0.1 (yield per recruit),
  calculated from the data frame in
  `Assessment@forecast[["per_recruit"]]`.

- `"Fmax"`: Fishing mortality relative to F_max (maximizing yield per
  recruit), calculated from the data frame in
  `Assessment@forecast[["per_recruit"]]`.

- `"FSPR"`: Fishing mortality relative to F_SPR% (the F that produces
  the spawning potential ratio specified in `"SPR_targ"`, calculated
  from data frame in `Assessment@forecast[["per_recruit"]]`.

- `"abs"`: Fishing mortality is independent of any model output and is
  explicitly specified in `relF`.

## Author

Q. Huynh

## Examples

``` r
# This is an MP with a 40-10 harvest control rule (using FMSY)
DD_40_10 <- make_MP(DD_TMB, HCR_segment, OCP_type = "SSB_SSB0", OCP = c(0.1, 0.4), relF = c(0, 1)) 
#' 
# This is an MP with a 40-10 harvest control rule with a maximum F of 0.1
DD_40_10 <- make_MP(DD_TMB, HCR_segment, OCP_type = "SSB_SSB0", 
                    Ftarget_type = "abs", OCP = c(0.1, 0.4), relF = c(0, 0.1)) 
```
