% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/RCM_data.R, R/define_S4_generic_RCM.R
\name{check_RCMdata}
\alias{check_RCMdata}
\alias{RCM}
\alias{RCM,OM,list-method}
\alias{RCM,OM,Data-method}
\alias{RCM,OM,RCMdata-method}
\alias{RCM,list,RCMdata-method}
\title{Rapid Conditioning Model (RCM)}
\usage{
check_RCMdata(RCMdata, OM, condition = "catch", silent = FALSE)

RCM(OM, data, ...)

\S4method{RCM}{OM,RCMdata}(
  OM,
  data,
  condition = "catch",
  selectivity = "logistic",
  s_selectivity = NULL,
  LWT = list(),
  comp_like = c("multinomial", "lognormal", "mvlogistic", "dirmult1", "dirmult2"),
  prior = list(),
  max_F = 3,
  cores = 1L,
  integrate = FALSE,
  mean_fit = FALSE,
  drop_nonconv = FALSE,
  drop_highF = FALSE,
  control = list(iter.max = 2e+05, eval.max = 4e+05),
  start = list(),
  map = list(),
  silent = FALSE,
  ...
)

\S4method{RCM}{OM,list}(
  OM,
  data,
  condition = "catch",
  selectivity = "logistic",
  s_selectivity = NULL,
  LWT = list(),
  comp_like = c("multinomial", "lognormal", "mvlogistic", "dirmult1", "dirmult2"),
  ESS = c(30, 30),
  prior = list(),
  max_F = 3,
  cores = 1L,
  integrate = FALSE,
  mean_fit = FALSE,
  drop_nonconv = FALSE,
  drop_highF = FALSE,
  control = list(iter.max = 2e+05, eval.max = 4e+05),
  start = list(),
  map = list(),
  silent = FALSE,
  ...
)

\S4method{RCM}{OM,Data}(
  OM,
  data,
  condition = "catch",
  selectivity = "logistic",
  s_selectivity = NULL,
  LWT = list(),
  comp_like = c("multinomial", "lognormal", "mvlogistic", "dirmult1", "dirmult2"),
  ESS = c(30, 30),
  prior = list(),
  max_F = 3,
  cores = 1L,
  integrate = FALSE,
  mean_fit = FALSE,
  drop_nonconv = FALSE,
  drop_highF = FALSE,
  control = list(iter.max = 2e+05, eval.max = 4e+05),
  start = list(),
  map = list(),
  silent = FALSE,
  ...
)

\S4method{RCM}{list,RCMdata}(
  OM,
  data,
  condition = "catch",
  selectivity = "logistic",
  s_selectivity = NULL,
  LWT = list(),
  comp_like = c("multinomial", "lognormal", "mvlogistic", "dirmult1", "dirmult2"),
  prior = list(),
  max_F = 3,
  integrate = FALSE,
  control = list(iter.max = 2e+05, eval.max = 4e+05),
  start = list(),
  map = list(),
  silent = FALSE,
  ...
)
}
\arguments{
\item{RCMdata}{An \linkS4class{RCMdata} object.}

\item{OM}{An object of class \link[MSEtool:OM-class]{MSEtool::OM} that specifies natural mortality (M), growth (Linf, K, t0, a, b), stock-recruitment relationship,
steepness, maturity parameters (L50 and L50_95), and standard deviation of recruitment variability (Perr).
Alternatively, provide a named list of biological inputs, see "StockPars" section below.}

\item{condition}{String to indicate whether the RCM is conditioned on "catch" (where F are estimated parameters), "catch2" (where F is solved internally using Newton's method),
or "effort" (F is proportional to an index series in \code{data@Ehist}. Can be fleet specific, in which case some combination of "catch" and "effort" are permissible.}

\item{silent}{Logical to indicate whether informative messages will be reported to console.}

\item{data}{Data inputs formatted in a \linkS4class{RCMdata} (preferred) or \link[MSEtool:Data-class]{MSEtool::Data} object.
Use of a list is deprecated. See Data section below.}

\item{...}{Other arguments to pass in for starting values of parameters and fixing parameters. See details.}

\item{selectivity}{A character vector of length nfleet to indicate \code{"logistic_length"}, \code{"dome_length"}, \code{"logistic_age"}, \code{"dome_age"}, or \code{"free"} selectivity for each fleet in \code{Chist}.
If there is time-varying selectivity, this is a character vector of length nsel_block (see Data section below). "free" indicates independent selectivity parameters for each age,
and additional modifications for fixing selectivity parameters will likely be needed. See Additional arguments section.}

\item{s_selectivity}{A vector of length nsurvey to indicate the selectivity of the corresponding columns in \code{data$Index}. Use \code{"B"} for
total biomass, or \code{"SSB"} for spawning biomass (by default, "B" is used). Use numbers if the survey selectivity follows a fleet (corresponding to the columns in data$Chist, e.g., 1 = first fleet/column and so on).
If the survey selectivity is otherwise independent of anything else in the model, use \code{"logistic_length"}, \code{"dome_length"}, \code{"logistic_age"}, \code{"dome_age"}, or \code{"free"} to specify the functional form of selectivity, and
see Additional arguments section for setup of survey selectivity parameters and Articles section for more information.}

\item{LWT}{A named list of likelihood weights for the RCM. See below.}

\item{comp_like}{A string indicating the statistical distribution for the composition data, either \code{"multinomial"} (default), \code{"lognormal"}, \code{"mvlogistic"} (multivariate logistic),
\code{"dirmult1"} (Dirichlet multinomial, linear version), or \code{"dirmult2"} (saturating version; see Thorson et al. 2017).}

\item{prior}{A named list for the parameters of any priors to be added to the model. See below.}

\item{max_F}{The maximum F for any fleet in the scoping model (higher F's in the model are penalized in the objective function). This argument will also update \code{OM@maxF}. See also \code{drop_highF}.}

\item{cores}{Integer for the number of CPU cores (set greater than 1 for parallel processing).}

\item{integrate}{Logical, whether to treat recruitment deviations as penalized parameters in the likelihood (FALSE) or random effects to be marginalized out of the likelihood (TRUE).}

\item{mean_fit}{Logical, whether to run an additional with mean values of life history parameters from the OM.}

\item{drop_nonconv}{Logical, whether to drop non-converged fits of the RCM, including fits where F = NA.}

\item{drop_highF}{Logical, whether to drop fits of the RCM where F = \code{max_F}.}

\item{control}{A named list of arguments (e.g, max. iterations, etc.) for optimization, to be passed to the control argument of \code{\link[stats:nlminb]{stats::nlminb()}}.}

\item{start}{A list of starting values for the TMB model. See details.}

\item{map}{A list of \code{map} argument to TMB models to override defaults. See \link[TMB:MakeADFun]{MakeADFun} and details.}

\item{ESS}{A vector of length two. A shortcut method to setting the maximum multinomial sample size of the age and length compositions.
Not used when data are provided in a \linkS4class{RCMdata} object.}
}
\value{
An object of class \linkS4class{RCModel} (see link for description of output).

\code{check_RCMdata} returns a list of updated RCMdata object, OM, and StockPars and FleetPars from the Hist object generated
from the OM.
}
\description{
Intended for conditioning operating models for MSEtool. For data-limited stocks, this function can generate a range of potential depletion scenarios inferred from sparse data.
From a historical time series of total catch or effort, and potentially age/length compositions and multiple indices of abundance, the RCM returns a range of values for depletion, selectivity,
unfished recruitment (R0), historical fishing effort, and recruitment deviations for the operating model. This is done by sampling life history parameters
provided by the user and fitting a statistical catch-at-age model (with the predicted catch equal to the observed catch).
Alternatively one can do a single model fit and sample the covariance matrix to generate an operating model with uncertainty based on the model fit.
Either a full catch (conditioned on catch) or effort (conditioned on effort) time series is needed but missing data (as NAs) are allowed for all other data types.
\code{check_RCMdata} evaluates whether the inputs in the S4 RCMdata object are correctly formatted.
}
\details{
Fleet selectivity is fixed to values sampled from \code{OM} if no age or length compositions are provided.

Survey selectivity is estimable only if \code{IAA} or \code{IAL} is provided. Otherwise, the selectivity should
be mirrored to a fleet (vulnerable biomass selectivity) or indexed to total or spawning biomass (see \code{s_selectivity}).

Parameters that were used in the fitting model are placed in the \code{RCM@OM@cpars} list.

If the operating model \code{OM} uses time-varying growth or M, then those trends will be used in the RCM as well.
Non-stationary productivity creates ambiguity in the calculation and interpretation of depletion and MSY reference points.

The easiest way to turn off time-varying growth/M is by setting: \code{OM@Msd <- OM@Linfsd <- OM@Ksd <- c(0, 0)}.

To play with alternative fits by excluding indices, for example, or other optional data, set the corresponding likelihood weight to zero. The model will still generate the inferred
index but the data won't enter the likelihood. See section on likelihood weights.
}
\section{Online Documentation}{

Several articles are available for RCM:

\itemize{
\item \href{https://openmse.com/tutorial-rcm/}{General overview of approach}
\item \href{https://openmse.com/tutorial-rcm-eq/}{Mathematical description}
\item \href{https://openmse.com/tutorial-rcm-select/}{Setup of selectivity settings and index catchability} (useful for more data-rich cases)
\item \href{https://openmse.com/features-assessment-models/5-priors/}{Description of priors}
}
}

\section{Priors}{

The following priors can be added as a named list, e.g., \verb{prior = list(M = c(0.25, 0.15), h = c(0.7, 0.1)}.
For each parameter below, provide a vector of values as described:

\describe{
\item{\code{R0}}{A vector of length 3. The first value indicates the distribution of the prior: \code{1} for lognormal, \code{2} for uniform
on \code{log(R0)}, \code{3} for uniform on R0. If lognormal, the second and third values are the prior mean (in normal space) and SD (in log space).
Otherwise, the second and third values are the lower and upper bounds of the uniform distribution (values in normal space).}
\item{\code{h}}{A vector of length 2 for the prior mean and SD, both in normal space. Beverton-Holt steepness uses a beta distribution,
while Ricker steepness uses a normal distribution.}
\item{\code{M}}{A vector of length 2 for the prior mean (in normal space) and SD (in log space). Lognormal prior.}
\item{\code{q}}{A matrix for nsurvey rows and 2 columns. The first column is the prior mean (in normal space) and the second column
for the SD (in log space). Use \code{NA} in rows corresponding to indices without priors.}
}
See online documentation for more details.
}

\section{Data}{

One of indices, age compositions, or length compositions should be provided in addition to the historical catch or effort. Not all arguments
are needed to run the model (some have defaults, while others are ignored if not applicable depending on the data provided).

The \code{data} variable can be an object of class \linkS4class{RCMdata}. See help file for description of inputs.

Alternatively, the \code{data} input can be a \link[MSEtool:Data-class]{MSEtool::Data} S4 object which will retrieve data from the following slots:

\describe{
\item{\code{Data@Cat}}{catch series (single fleet with the Data S4 object)}
\item{\code{Data@Effort}}{effort series}
\item{\code{Data@CAA}}{fishery age composition}
\item{\code{Data@CAL}, \code{Data@CAL_mids}}{fishery length composition and corresponding length bins}
\item{\code{Data@Ind}, \code{Data@SpInd}, \code{Data@VInd}, \code{Data@AddInd}}{indices of abundance}
\item{\code{Data@CV_Ind}, \code{Data@CV_SpInd}, \code{Data@CV_VInd}, \code{Data@CV_AddInd}}{annual coefficients of variation for the corresponding indices
of abundance. CVs will be converted to lognormal standard deviations.}
\item{\code{Data@ML}}{fishery mean lengths}
\item{\code{Data@AddIndV}, \code{Data@AddIndType}, \code{Data@AddIunits}}{Additional information for indices in \code{Data@AddInd}:
selectivity and units (i.e., biomass or abundance).}
}

There is no slot in the Data S4 object for the equilibrium catch/effort. These can be passed directly in the function call, i.e., \code{RCM(OM, Data, C_eq = C_eq, ...)}.
}

\section{Data list (deprecated)}{

Use of a list is deprecated. For backwards compatibility, here is the list of supported entries:

\describe{
\item{\code{Chist}}{A vector of historical catch, should be of length OM@nyears. If there are multiple fleets: a matrix of \code{OM@nyears} rows and \code{nfleet} columns.
Ideally, the first year of the catch series represents unfished conditions (see also \code{C_eq}).}
\item{\code{C_sd}}{A vector or matrix of standard deviations (lognormal distribution) for the catches in \code{Chist}.
If not provided, the default is 0.01. Only used if \code{condition = "catch"}.}
\item{\code{Ehist}}{A vector of historical effort, should be of length \code{OM@nyears} (see also \code{E_eq}).}
\item{\code{Index}}{A vector of values of an index (of length \code{OM@nyears}). If there are multiple indices: a matrix of historical indices of abundances, with rows
indexing years and columns indexing the index.}
\item{\code{I_sd}}{A vector or matrix of standard deviations (lognormal distribution) for the indices corresponding to the entries in \code{Index}.
If not provided, this function will use values from \code{OM@Iobs}.}
\item{\code{I_type}}{Obsolete as of version 2.0. See \code{s_selectivity} argument.}
\item{\code{CAA}}{Fishery age composition matrix with \code{nyears} rows and \code{OM@maxage+1} columns. If multiple fleets: an array with dimension:
\verb{nyears, OM@maxage, and nfleet}.}
\item{\code{CAL}}{Fishery length composition matrix with nyears rows and columns indexing the length bin. If multiple fleets:
an array with dimension: \verb{nyears, length bins, and nfleet}.}
\item{\code{MS}}{A vector of fishery mean size (MS, either mean length or mean weight) observations (length \code{OM@nyears}),
or if multiple fleets: matrix of dimension: \verb{nyears, nfleet}.
Generally, mean lengths should not be used if \code{CAL} is also provided, unless mean length and length comps are independently sampled.}
\item{\code{MS_type}}{A character (either \code{"length"} (default) or \code{"weight"}) to denote the type of mean size data.}
\item{\code{MS_cv}}{The coefficient of variation of the observed mean size. If there are multiple fleets, a vector of length \code{nfleet}.
Default is 0.2.}
\item{\code{s_CAA}}{Survey age composition data, an array of dimension \verb{nyears, maxage+1, nsurvey}.}
\item{\code{s_CAL}}{Survey length composition data, an array of dimension \verb{nyears, length(length_bin), nsurvey}.}
\item{\code{length_bin}}{A vector for the midpoints of the length bins for \code{CAL} and \code{s_CAL}. All bin widths should be equal in size.}
\item{\code{C_eq}}{A numeric vector of length \code{nfleet} for the equilibrium catch for each fleet in \code{Chist} prior to the first year of the operating model.
Zero (default) implies unfished conditions in year one. Otherwise, this is used to estimate depletion in the first year of the data. Alternatively,
if one has a full CAA matrix, one could instead estimate "artificial" rec devs to generate the initial numbers-at-age (and hence initial
depletion) in the first year of the model (see additional arguments).}
\item{\code{C_eq_sd}}{A vector of standard deviations (lognormal distribution) for the equilibrium catches in \code{C_eq}.
If not provided, the default is 0.01. Only used if \code{condition = "catch"}.}
\item{\code{E_eq}}{The equilibrium effort for each fleet in \code{Ehist} prior to the first year of the operating model.
Zero (default) implies unfished conditions in year one. Otherwise, this is used to estimate depletion in the first year of the data.}
\item{\code{abs_I}}{Optional, an integer vector to indicate which indices are in absolute magnitude. Use 1 to set \code{q = 1},
otherwise use 0 to estimate q.}
\item{\code{I_units}}{Optional, an integer vector to indicate whether indices are biomass based (1) or abundance-based (0).
By default, all are biomass-based.}
\item{\code{age_error}}{Optional, a square matrix of maxage + 1 rows and columns to specify ageing error. The aa-th column
assigns a proportion of the true age in the a-th row to observed age. Thus, all rows should sum to 1.
Default is an identity matrix (no ageing error).}
\item{\code{sel_block}}{Optional, for time-varying fleet selectivity (in time blocks), a integer matrix of \code{nyears} rows and \code{nfleet} columns
to assigns a selectivity function to a fleet for certain years.}
}
}

\section{StockPars}{

When an operating model is provided, the RCM function will generally fit to each simulation of biological parameters.

Alternatively for a single fit to data independent of any operating model, provide a named list containing the following
(naming conventions follow internal operating model variables):
\itemize{
\item \code{SRrel} Integer, stock-recruit function (1 = Beverton-Holt, 2 = Ricker, 3 = Mesnil-Rochet hockey stick)
\item \code{R0} Numeric, starting value for unfished recruitment parameter
\item \code{M_ageArray} Matrix \verb{[maxage+1, nyears]} for natural mortality
\item \code{Len_age} Matrix \verb{[maxage+1, nyears + 1]} for length at age
\item \code{Linf} Numeric. Asymptotic length. Only used for the upper bound for the size of full selectivity (if selectivity functions are length-based)
\item \code{LatASD} Matrix \verb{[maxage+1, nyears + 1]} for the standard deviation in length at age
\item \code{Wt_age} Matrix \verb{[maxage+1, nyears + 1]} for stock weight at age
\item \code{Mat_age} Matrix \verb{[maxage+1, nyears + 1]} for maturity at age
\item \code{Fec_Age} Matrix \verb{[maxage+1, nyears + 1]} for fecundity at age. Frequently the product of maturity and weight at age
\item \code{ageMarray} Numeric, age of 50 percent maturity. Used to average the initial years for the unfished replacement line of the
stock recruit relationship and steepness/R0. Irrelevant if fecundity and natural mortality are not time-varying (set to 1).
\item \code{spawn_time_frac} Numeric, fraction of the year when spawning occurs
\item \code{hs} Numeric, steepness of the stock recruit relationship
\item \code{procsd} Numeric, lognormal recruitment deviation standard deviation
}
}

\section{Additional arguments}{

For \code{RCM}, additional arguments can be passed to the model via \code{...}:

\describe{
\item{\code{plusgroup}}{Logical for whether the maximum age is a plusgroup or not. By default, TRUE.}
\item{\code{fix_dome}}{Logical for whether the dome selectivity parameter for fleets is fixed. Used primarily for backwards compatibility,
this is overridden by the \code{map} argument.}
\item{\code{resample}}{Logical, whether the OM conditioning parameters (recruitment, fishing mortality, SSB, selectivity, etc.) are obtained by sampling the Hessian matrix from
a single model fit. By default FALSE. This feature requires identical biological parameters among simulations.}
\item{\code{pbc_recdev}}{Vector of length nyears. Proportion of the bias correction to apply annually to the recruitment deviations (if estimated).
The bias correction from logspace to normal space is \code{exp(log_rec_dev[y] - 0.5 * pbc_recdev[y] * sigmaR^2)}. Default proportion is 1.}
\item{\code{pbc_earlyrecdev}}{Vector of length maxage. Proportion of the bias correction to apply to the abundance deviations in the first year of the model (if estimated).
The bias correction from logspace to normal space is \code{exp(log_early_rec_dev[a] - 0.5 * pbc_recdev[a] * sigmaR^2)}. Default proportion is 1.}
}
}

\section{start}{

Starting values can be specified in a named list for the following:

\describe{
\item{\code{vul_par}}{A matrix of 3 rows and nfleet columns for starting values for fleet selectivity. The three rows correspond
to LFS (length of full selectivity), L5 (length of 5 percent selectivity), and Vmaxlen (selectivity at length Linf). By default,
the starting values are values from the OM object. If any \code{selectivity = "free"}, then this matrix needs to be of \code{maxage+1} rows where
the row specifies the selectivity at age. See Articles section.}
\item{\code{ivul_par}}{A matrix of 3 rows and nsurvey columns for starting values for fleet selectivity. Same setup as \code{vul_par}. Values in the column are ignored
if \code{s_selectivity} is mapped to a fishing fleet (add NA placeholders in that case).
If any \code{s_selectivity = "free"}, then this matrix needs to be of \code{maxage+1} rows where
the row specifies the selectivity at age.}
\item{\code{log_rec_dev}}{A numeric vector of length \code{nyears} for the starting values of the log-recruitment deviations.}
\item{\code{log_early_rec_dev}}{A numeric vector of length \code{OM@maxage} for the starting values of the recruitment deviations controlling the abundance-at-age in the first year of the model.}
\item{\code{q}}{A numeric vector of length nsurvey for index catchability. See \href{https://openmse.com/tutorial-rcm-select/}{online article} for more information.}
}
}

\section{map}{

Parameters can be fixed with the map argument (also a named list, corresponding to the start list). Each
vector or matrix in the map argument will be the same dimension as in the start entry. If an entry is \code{NA}, the corresponding parameter is fixed in the model to the starting
value. Otherwise, an integer for each independent parameter, i.e., shared or mirrored parameters get the same integer entry.

\describe{
\item{\code{vul_par}}{An integer matrix of the same dimension as \code{start$vul_par}. By default, selectivity is fixed if there are no age or length composition for that fleet
or survey, otherwise estimated. Unused cells in the \code{start$vul_par} matrix should be given NA in the map matrix.}
\item{\code{ivul_par}}{The map argument for the survey selectivity parameters (same dimension as \code{start$ivul_par}). Placeholder parameters should have a map value of NA.}
\item{\code{log_early_rec_dev}}{A vector of length \code{OM@maxage} that indexes which recruitment deviates for the cohorts in the first year of the model are fixed (using NA) or estimated (a separate integer).
By default, no deviates are estimated (all are NA).}
\item{\code{log_rec_dev}}{A vector of length \code{OM@nyears} that indexes which recruitment deviates are fixed (using NA) or estimated (a separate integer).
By default, all these deviates are estimated.}
\item{\code{q}}{A vector of length \code{nsurvey} for index catchability. q should be an estimated parameter when sharing across surveys (perhaps with differing selectivity). Otherwise, it is solved analytically
where individual parameters are independent of other indices. Use \code{RCMdata@abs_I} for fixing the catchability to 1. See \href{https://openmse.com/tutorial-rcm-select/}{online article} for more information.}
}
}

\section{Likelihood weights}{

\code{LWT} is an optional named list containing the likelihood weights (values >= 0) with the possible options:
\itemize{
\item \verb{Chist, CAA, CAL, MS, C_eq}: A vector of length nfleet for each.
\item \verb{Index, IAA, IAL}: A vector of length nsurvey for each.
}

By default, all likelihood weights are equal to one if not specified by the user.

Annual multinomial sample sizes for the age and length comps can now be provided directly in the
\linkS4class{RCMdata} object. For a list or \link[MSEtool:Data-class]{MSEtool::Data} object, use the \code{ESS} argument.
}

\examples{
\donttest{ 
# An example that conditions a Pacific cod operating model. There are 48 simulations, 
# where values of natural mortality and steepness are sampled from distributions. 
# The model is fitted with priors on the index catchability. Maturity and selectivity 
# are knife-edge at the age of 2 years. See online tutorial for more information.

data(pcod) 
mat_ogive <- pcod$OM@cpars$Mat_age[1, , 1]
out <- RCM(OM = pcod$OM, data = pcod$data, 
           condition = "catch", mean_fit = TRUE,
           selectivity = "free", s_selectivity = rep("SSB", ncol(pcod$data@Index)),
           start = list(vul_par = matrix(mat_ogive, length(mat_ogive), 1)),
           map = list(vul_par = matrix(NA, length(mat_ogive), 1),
                      log_early_rec_dev = rep(1, pcod$OM@maxage)),
           prior = pcod$prior)
plot(out, s_name = colnames(pcod$data@Index))

# Alternative OM with age-3 maturity and selectivity instead.
out_age3 <- local({
  pcod$OM@cpars$Mat_age[, 2, ] <- 0
  mat_ogive_age3 <- pcod$OM@cpars$Mat_age[1, , 1]
  RCM(OM = pcod$OM, data = pcod$data, 
      condition = "catch", mean_fit = TRUE,
      selectivity = "free", s_selectivity = rep("SSB", ncol(pcod$data@Index)),
      start = list(vul_par = matrix(mat_ogive_age3, length(mat_ogive_age3), 1)),
      map = list(vul_par = matrix(NA, length(mat_ogive_age3), 1),   
                 log_early_rec_dev = rep(1, pcod$OM@maxage)),
      prior = pcod$prior)
})
  
compare_RCM(out, out_age3, scenario = list(names = c("Age-2 maturity", "Age-3 maturity")),
            s_name = colnames(pcod$data@Index))
             
Hist <- runMSE(out@OM, Hist = TRUE)            
} 

}
\references{
Thorson et al. 2017. Model-based estimates of effective sample size in stock assessment models using the Dirichlet-multinomial distribution.
Fish. Res. 192:84-93. \doi{10.1016/j.fishres.2016.06.005}
}
\seealso{
\link{plot.RCModel} \linkS4class{RCModel} \link{compare_RCM} \link{pcod} \link{RCM2MOM} \link{posterior}
}
\author{
Q. Huynh
}
