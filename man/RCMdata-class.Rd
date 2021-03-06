% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/define_S4_class_assess_RCMdata.R
\docType{class}
\name{RCMdata-class}
\alias{RCMdata-class}
\alias{RCMdata}
\title{Class-\code{RCMdata}}
\description{
An S4 class for the data inputs into \link{RCM}.
}
\section{Slots}{

\describe{
\item{\code{Chist}}{A vector of historical catch, should be of length OM@nyears. If there are multiple fleets: a matrix of OM@nyears rows and nfleet columns.
Ideally, the first year of the catch series represents unfished conditions (see also \code{C_eq}).}

\item{\code{C_sd}}{A vector or matrix of lognormal distribution standard deviations (by year and fleet) for the catches in \code{Chist}.
If not provided, the default is 0.01. Not used if \code{RCM(condition = "catch2")}.}

\item{\code{Ehist}}{A vector of historical effort, should be of length OM@nyears. If there are multiple fleets: a matrix of OM@nyears rows and nfleet columns. 
See also \code{E_eq}).}

\item{\code{CAA}}{Fishery age composition matrix with nyears rows and OM@maxage+1 columns. If multiple fleets: an array with dimension: nyears, OM@maxage, and nfleets.
Raw numbers will be converted to annual proportions (see CAA_ESS for sample sizes).}

\item{\code{CAA_ESS}}{Annual sample size (for the multinomial distribution) of the fishery age comps. 
A vector of length OM@nyears. If there are multiple fleets: a matrix of OM@nyears rows and nfleet columns.}

\item{\code{CAL}}{Fishery length composition matrix with nyears rows and columns indexing the length bin. If multiple fleets: an array with dimension: nyears,
n_bins, and nfleets. Raw numbers will be converted to annual proportions (see CAL_ESS for sample sizes).}

\item{\code{CAL_ESS}}{Annual sample size (for the multinomial distribution) of the fishery length comps. 
A vector of length OM@nyears. If there are multiple fleets: a matrix of OM@nyears rows and nfleet columns.}

\item{\code{length_bin}}{- A vector (length n_bin) for the midpoints of the length bins for \code{CAL} and \code{IAL}, as well as the population model, if all bin widths are equal in size. 
If length bins are unequal in bin width, then provide a vector of the boundaries of the length bins (vector of length n_bin + 1).}

\item{\code{MS}}{A vector of fishery mean size (MS, either mean length or mean weight) observations (length OM@nyears), or if multiple fleets: matrix of dimension: nyears and nfleets.
Generally, mean lengths should not be used if \code{CAL} is also provided, unless mean length and length comps are independently sampled.}

\item{\code{MS_type}}{A character (either \code{"length"} (default) or \code{"weight"}) to denote the type of mean size data.}

\item{\code{MS_cv}}{The coefficient of variation of the observed mean size. If there are multiple fleets, a vector of length nfleet.
Default is 0.2.}

\item{\code{Index}}{A vector of values of an index (of length OM@nyears). If there are multiple surveys: a matrix of historical indices of abundances, with rows
indexing years and columns indexing surveys.}

\item{\code{I_sd}}{A vector or matrix of standard deviations (lognormal distribution) for the indices corresponding to the entries in \code{Index}.
If not provided, this function will use values from \code{OM@Iobs}.}

\item{\code{IAA}}{Index age composition data, an array of dimension nyears, maxage+1, nsurvey.
Raw numbers will be converted to annual proportions (see IAA_ESS for sample sizes).}

\item{\code{IAA_ESS}}{Annual sample size (for the multinomial distribution) of the index age comps. 
A vector of length OM@nyears. If there are multiple indices: a matrix of OM@nyears rows and nsurvey columns.}

\item{\code{IAL}}{Index length composition data, an array of dimension nyears, n_bin, nsurvey. 
Raw numbers will be converted to annual proportions (see IAL_ESS for sample sizes).}

\item{\code{IAL_ESS}}{Annual sample size (for the multinomial distribution) of the index length comps. 
A vector of length OM@nyears. If there are multiple indices: a matrix of OM@nyears rows and nsurvey columns.}

\item{\code{C_eq}}{A numeric vector of length nfleet for the equilibrium catch for each fleet in \code{Chist} prior to the first year of the operating model.
Zero (default) implies unfished conditions in year one. Otherwise, this is used to estimate depletion in the first year of the data. Alternatively,
if one has a full CAA matrix, one could instead estimate "artificial" rec devs to generate the initial numbers-at-age (and hence initial depletion) 
in the first year of the model (see additional arguments in \link{RCM}).}

\item{\code{C_eq_sd}}{- A vector of standard deviations (lognormal distribution) for the equilibrium catches in \code{C_eq}.
If not provided, the default is 0.01. Only used if \code{RCM(condition = "catch")}.}

\item{\code{E_eq}}{The equilibrium effort for each fleet in \code{Ehist} prior to the first year of the operating model.
Zero (default) implies unfished conditions in year one. Otherwise, this is used to estimate depletion in the first year of the data.}

\item{\code{abs_I}}{An integer vector to indicate which indices are in absolute magnitude. Use 1 to set q = 1, otherwise use 0 (default) to estimate q.}

\item{\code{I_units}}{An integer vector to indicate whether indices are biomass based (1) or abundance-based (0). By default, all are biomass-based.}

\item{\code{age_error}}{A square matrix of maxage + 1 rows and columns to specify ageing error. The aa-th column assigns a proportion of the true age in the
a-th row to observed age. Thus, all rows should sum to 1. Default is an identity matrix (no ageing error).}

\item{\code{sel_block}}{For time-varying fleet selectivity (in time blocks), a integer matrix of nyears rows and nfleet columns to assigns a selectivity function 
to a fleet for certain years. By default, constant selectivity for each individual fleet.
See the \href{https://openmse.com/tutorial-rcm-select/}{selectivity} article for more details.}

\item{\code{Misc}}{A list of miscellaneous inputs. Used internally.}
}}

\seealso{
\link{RCM}
}
\author{
Q. Huynh
}
