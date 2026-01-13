# The rapid conditioning model as an assessment function

In beta testing. A function that uses RCM as an assessment function for
use in MPs. More function arguments will be added to tinker with model
settings and data inputs.

## Usage

``` r
RCM_assess(
  x = 1,
  Data,
  AddInd = "B",
  SR = c("BH", "Ricker"),
  selectivity = c("logistic", "dome"),
  CAA_multiplier = 50,
  prior = list(),
  LWT = list(),
  StockPars = "Data",
  ...
)
```

## Arguments

- x:

  A position in the Data object (by default, equal to one for
  assessments).

- Data:

  An object of class Data

- AddInd:

  A vector of integers or character strings indicating the indices to be
  used in the model. Integers assign the index to the corresponding
  index in Data@AddInd, "B" (or 0) represents total biomass in Data@Ind,
  "VB" represents vulnerable biomass in Data@VInd, and "SSB" represents
  spawning stock biomass in Data@SpInd. Vulnerability to the survey is
  fixed in the model.

- SR:

  Stock-recruit function (either `"BH"` for Beverton-Holt or
  `"Ricker"`).

- selectivity:

  Whether to model "logistic" or "dome" selectivity for the fishery.

- CAA_multiplier:

  Numeric for data weighting of catch-at-age matrix. If greater than 1,
  then this is the maximum multinomial sample size in any year. If less
  than one, then the multinomial sample size is this fraction of the
  sample size.

- prior:

  A named list for the parameters of any priors to be added to the
  model. See documentation in
  [SCA](https://samtool.openmse.com/reference/SCA.md).

- LWT:

  A named list (Index, CAA, Catch) of likelihood weights for the data
  components. For the index, a vector of length survey. For CAL and
  Catch, a single value.

- StockPars:

  Either a string ("Data" or "OM") to indicate whether to grab
  biological parameters from the Data object, or operating model.
  Alternatively, a named list to provide custom parameters for the
  assessment.

- ...:

  Additional arguments (to be added).

## Details

*Data* object: the function currently uses catch, CAA, and indices of
abundance in the corresponding slots in the Data object.

*StockPars* input: biological parameters can be used from (1) Data
object, (2) operating model, or (3) provided directly in the `StockPars`
argument.

Options 2 and 3 allow for time-varying growth, maturity, and natural
mortality. Natural mortality can also be age-varying.

`StockPars` can be a named list of parameters used to provide inputs to
the assessment model:

- `Wt_age` - annual weight at age, array `[sim, ages, year]`

- `Mat_age` - annual maturity at age, array `[sim, ages, year]`

- `hs` - Stock-recruit steepness, vector of length `[sim]`

- `M_ageArray` - annual natural mortality, array `[sim, ages, year]`

## Examples

``` r
 
r <- RCM_assess(Data = SimulatedData)
myMP <- make_MP(RCM_assess, HCR_MSY)
myMP(x = 1, Data = SimulatedData)
#> TAC (median) 
#>     23620.79 
```
