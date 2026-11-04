
<!-- README.md is generated from README.Rmd. Please edit that file -->

# SAMtool <img src='man/figures/logo.png' align="right" height="139" />

<!-- badges: start -->

[![CRAN
status](https://www.r-pkg.org/badges/version/SAMtool)](https://CRAN.R-project.org/package=SAMtool)
[![R-CMD-check](https://github.com/Blue-Matter/SAMtool/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/Blue-Matter/SAMtool/actions/workflows/R-CMD-check.yaml)
<!-- badges: end -->

`SAMtool` is the Stock Assessment Methods toolkit, a collection of stock
assessment models and assessment-based management procedures for
data-moderate and data-rich fisheries. It also provides the Rapid
Conditioning Model (RCM) for conditioning operating models on fishery
data.

`SAMtool` is one of the three core packages of the
[openMSE](https://openMSE.com) collection, alongside
[MSEtool](https://msetool.openmse.com), which provides operating models
and closed-loop simulation, and [DLMtool](https://dlmtool.openmse.com),
which provides data-limited management procedures. Installing `openMSE`
installs all three, and `library(openMSE)` attaches them.

## Key Features

- **Operating model conditioning**: the Rapid Conditioning Model
  (`RCM()`) fits a statistical catch-at-age model to catch or effort,
  indices of abundance, and age and length compositions, and returns an
  operating model conditioned on the fit
- **Conditioning of MSEtool operating models**: `RCM2()` takes `stock`,
  `fleet`, and `data` objects, fits each stock (including hermaphroditic
  stocks) with fleet-specific selectivity and fishing mortality, and
  returns an `om` object ready for `Simulate()` and `runMSE()`
- **Assessment models**: statistical catch-at-age (`SCA()`),
  delay-difference (`DD_TMB()`, `cDD()`), surplus production (`SP()`),
  virtual population analysis (`VPA()`), and simple stock synthesis
  (`SSS()`) models, built with
  [TMB](https://CRAN.R-project.org/package=TMB)
- **Assessment-based management procedures**: `make_MP()` combines any
  assessment model with a harvest control rule, and many ready-made MPs
  are included, along with interim and projection-based MPs
- **Harvest control rules**: ramped, threshold, fixed-F, escapement, and
  MSY-based HCRs that can be combined with any assessment model
- **Diagnostics**: standardized reports, retrospective analyses,
  likelihood profiles, and diagnostics of assessment performance inside
  closed-loop simulation
- **Bayesian estimation**: MCMC sampling of RCM and assessment model
  posteriors with [tmbstan](https://CRAN.R-project.org/package=tmbstan)

## Installation

Install the released version from [CRAN](https://CRAN.R-project.org):

``` r
install.packages("SAMtool")
```

Or the development version from GitHub:

``` r
# install.packages("pak")
pak::pak("Blue-Matter/SAMtool")
```

## Documentation

- **[openMSE.com](https://openMSE.com)**: general information about the
  openMSE framework
- **[docs.openMSE.com](https://docs.openMSE.com)**: technical manual
  with worked examples covering operating model construction,
  simulation, and MSE analysis
- **[samtool.openmse.com](https://samtool.openmse.com)**: full function
  reference for `SAMtool`
