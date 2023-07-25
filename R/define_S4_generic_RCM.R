

#' @name RCM
#' @aliases RCM RCM,OM,list-method RCM,OM,Data-method
#' @title Rapid Conditioning Model (RCM)
#'
#' @description Intended for conditioning operating models for MSEtool. For data-limited stocks, this function can generate a range of potential depletion scenarios inferred from sparse data.
#' From a historical time series of total catch or effort, and potentially age/length compositions and multiple indices of abundance, the RCM returns a range of values for depletion, selectivity,
#' unfished recruitment (R0), historical fishing effort, and recruitment deviations for the operating model. This is done by sampling life history parameters
#' provided by the user and fitting a statistical catch-at-age model (with the predicted catch equal to the observed catch).
#' Alternatively one can do a single model fit and sample the covariance matrix to generate an operating model with uncertainty based on the model fit.
#' Either a full catch (conditioned on catch) or effort (conditioned on effort) time series is needed but missing data (as NAs) are allowed for all other data types.
#' \code{check_RCMdata} evaluates whether the inputs in the S4 RCMdata object are correctly formatted.
#' 
#' @param OM An object of class \linkS4class{OM} that specifies natural mortality (M), growth (Linf, K, t0, a, b), stock-recruitment relationship,
#' steepness, maturity parameters (L50 and L50_95), standard deviation of recruitment variability (Perr), as well as index uncertainty (Iobs).
#' @param data Data inputs formatted in a \linkS4class{RCMdata} (preferred) or \linkS4class{Data} object. 
#' Use of a list is deprecated. See Data section below.
#' @param condition String to indicate whether the RCM is conditioned on "catch" (where F are estimated parameters), "catch2" (where F is solved internally using Newton's method),
#' or "effort" (F is proportional to an index series in \code{data@Ehist}. Can be fleet specific.
#' @param selectivity A character vector of length nfleet to indicate \code{"logistic_length"}, \code{"dome_length"}, \code{"logistic_age"}, \code{"dome_age"}, or \code{"free"} selectivity for each fleet in \code{Chist}.
#' If there is time-varying selectivity, this is a character vector of length nsel_block (see Data section below). "free" indicates independent selectivity parameters for each age,
#' and additional modifications for fixing selectivity parameters will likely be needed. See Additional arguments section.
#' @param s_selectivity A vector of length nsurvey to indicate the selectivity of the corresponding columns in \code{data$Index}. Use \code{"B"} for
#' total biomass, or \code{"SSB"} for spawning biomass (by default, "B" is used). Use numbers if the survey selectivity follows a fleet (corresponding to the columns in data$Chist, e.g., 1 = first fleet/column and so on).
#' If the survey selectivity is otherwise independent of anything else in the model, use \code{"logistic_length"}, \code{"dome_length"}, \code{"logistic_age"}, \code{"dome_age"}, or \code{"free"} to specify the functional form of selectivity, and
#' see Additional arguments section for setup of survey selectivity parameters and Articles section for more information.
#' @param LWT A named list of likelihood weights for the RCM. See below.
#' @param comp_like A string indicating the statistical distribution for the composition data, either \code{"multinomial"} (default), \code{"lognormal"}, \code{"mvlogistic"} (multivariate logistic),
#' \code{"dirmult1"} (Dirichlet multinomial, linear version), or \code{"dirmult2"} (saturating version; see Thorson et al. 2017).
#' @param ESS A vector of length two. A shortcut method to setting the maximum multinomial sample size of the age and length compositions. 
#' Not used when data are provided in a \linkS4class{RCMdata} object.
#' @param prior A named list for the parameters of any priors to be added to the model. See below.
#' @param max_F The maximum F for any fleet in the scoping model (higher F's in the model are penalized in the objective function). This argument will also update \code{OM@maxF}. See also \code{drop_highF}.
#' @param cores Integer for the number of CPU cores (set greater than 1 for parallel processing).
#' @param integrate Logical, whether to treat recruitment deviations as penalized parameters in the likelihood (FALSE) or random effects to be marginalized out of the likelihood (TRUE).
#' @param mean_fit Logical, whether to run an additional with mean values of life history parameters from the OM.
#' @param drop_nonconv Logical, whether to drop non-converged fits of the RCM, including fits where F = NA.
#' @param drop_highF Logical, whether to drop fits of the RCM where F = \code{max_F}.
#' @param control A named list of arguments (e.g, max. iterations, etc.) for optimization, to be passed to the control argument of \code{\link[stats]{nlminb}}.
#' @param start A list of starting values for the TMB model. See details.
#' @param map A list of \code{map} argument to TMB models to override defaults. See \link[TMB]{MakeADFun} and details.
#' @param silent Logical to indicate whether informative messages will be reported to console.
#' @param ... Other arguments to pass in for starting values of parameters and fixing parameters. See details.
#'
#' @section Priors:
#' The following priors can be added as a named list, e.g., \code{prior = list(M = c(0.25, 0.15), h = c(0.7, 0.1)}. 
#' For each parameter below, provide a vector of values as described:
#' 
#' \itemize{
#' \item \code{R0} - A vector of length 3. The first value indicates the distribution of the prior: \code{1} for lognormal, \code{2} for uniform
#' on \code{log(R0)}, \code{3} for uniform on R0. If lognormal, the second and third values are the prior mean (in normal space) and SD (in log space).
#' Otherwise, the second and third values are the lower and upper bounds of the uniform distribution (values in normal space).
#' \item \code{h} - A vector of length 2 for the prior mean and SD, both in normal space. Beverton-Holt steepness uses a beta distribution, 
#' while Ricker steepness uses a normal distribution.
#' \item \code{M} - A vector of length 2 for the prior mean (in normal space) and SD (in log space). Lognormal prior.
#' \item \code{q} - A matrix for nsurvey rows and 2 columns. The first column is the prior mean (in normal space) and the second column 
#' for the SD (in log space). Use \code{NA} in rows corresponding to indices without priors.
#' }
#' See online documentation for more details.
#' 
#' @details
#' Fleet selectivity is fixed to values sampled from \code{OM} if no age or length compositions are provided.
#'
#' Survey selectivity is estimable only if \code{IAA} or \code{IAL} is provided. Otherwise, the selectivity should
#' be mirrored to a fleet (vulnerable biomass selectivity) or indexed to total or spawning biomass (see \code{s_selectivity}).
#'
#' Parameters that were used in the fitting model are placed in the \code{RCM@@OM@@cpars} list.
#' 
#' If the operating model \code{OM} uses time-varying growth or M, then those trends will be used in the RCM as well.
#' Non-stationary productivity creates ambiguity in the calculation and interpretation of depletion and MSY reference points.
#'
#' The easiest way to turn off time-varying growth/M is by setting: \code{OM@@Msd <- OM@@Linfsd <- OM@@Ksd <- c(0, 0)}.
#'
#' To play with alternative fits by excluding indices, for example, or other optional data, set the corresponding likelihood weight to zero. The model will still generate the inferred
#' index but the data won't enter the likelihood. See section on likelihood weights.
#'
#' @return An object of class \linkS4class{RCModel} (see link for description of output).
#' 
#' \code{check_RCMdata} returns a list of updated RCMdata object, OM, and StockPars and FleetPars from the Hist object generated
#' from the OM.
#'
#' @section Online Documentation:
#' Several articles are available for the RCM:
#'
#' \itemize{
#' \item \href{https://openmse.com/tutorial-rcm/}{General overview of approach}
#' \item \href{https://openmse.com/tutorial-rcm-eq/}{Mathematical description}
#' \item \href{https://openmse.com/tutorial-rcm-select/}{Setup of selectivity settings and index catchability} (useful for more data-rich cases)
#' \item \href{https://openmse.com/features-assessment-models/5-priors/}{Description of priors}
#' }
#'
#' @section Data:
#' One of indices, age compositions, or length compositions should be provided in addition to the historical catch or effort. Not all arguments
#' are needed to run the model (some have defaults, while others are ignored if not applicable depending on the data provided).
#'
#' The \code{data} variable can be an object of class \linkS4class{RCMdata}. See help file for description of inputs.
#'
#' Alternatively, the \code{data} input can be a \linkS4class{Data} S4 object which will retrieve data from the following slots:
#'
#' \itemize{
#' \item Data@@Cat - catch series (single fleet with the Data S4 object)
#' \item Data@@Effort - effort series
#' \item Data@@CAA - fishery age composition
#' \item Data@@CAL, Data@@CAL_mids - fishery length composition and corresponding length bins
#' \item Data@@Ind, Data@@SpInd, Data@@VInd, Data@@AddInd - indices of abundance
#' \item Data@@CV_Ind, Data@@CV_SpInd, Data@@CV_VInd, Data@@CV_AddInd - annual coefficients of variation for the corresponding indices of abundance. CVs will be converted to lognormal standard deviations.
#' \item Data@@ML - fishery mean lengths
#' \item Data@@AddIndV, Data@@AddIndType, Data@@AddIunits - Additional information for indices in Data@@AddInd: selectivity and units (i.e., biomass or abundance).
#' }
#'
#' There is no slot in the Data S4 object for the equilibrium catch/effort. These can be passed directly in the function call, i.e., \code{RCM(OM, Data, C_eq = C_eq, ...)}.
#'
#' Use of a list is deprecated. For backwards compatibility, here is the list of supported entries:
#' 
#' \itemize{
#' \item Chist - A vector of historical catch, should be of length OM@@nyears. If there are multiple fleets: a matrix of OM@@nyears rows and nfleet columns.
#' Ideally, the first year of the catch series represents unfished conditions (see also \code{C_eq}).
#' \item C_sd - A vector or matrix of standard deviations (lognormal distribution) for the catches in \code{Chist}.
#' If not provided, the default is 0.01. Only used if \code{condition = "catch"}.
#' \item Ehist - A vector of historical effort, should be of length OM@@nyears (see also \code{E_eq}).
#' \item Index - A vector of values of an index (of length OM@@nyears). If there are multiple indices: a matrix of historical indices of abundances, with rows
#' indexing years and columns indexing the index.
#' \item I_sd - A vector or matrix of standard deviations (lognormal distribution) for the indices corresponding to the entries in \code{Index}.
#' If not provided, this function will use values from \code{OM@@Iobs}.
#' \item I_type - Obsolete as of version 2.0. See \code{s_selectivity} argument.
#' \item CAA - Fishery age composition matrix with nyears rows and OM@@maxage+1 columns. If multiple fleets: an array with dimension: nyears, OM@@maxage, and nfleets.
#' \item CAL - Fishery length composition matrix with nyears rows and columns indexing the length bin. If multiple fleets: an array with dimension: nyears,
#' length bins, and nfleets.
#' \item MS - A vector of fishery mean size (MS, either mean length or mean weight) observations (length OM@@nyears), or if multiple fleets: matrix of dimension: nyears and nfleets.
#' Generally, mean lengths should not be used if \code{CAL} is also provided, unless mean length and length comps are independently sampled.
#' \item MS_type - A character (either \code{"length"} (default) or \code{"weight"}) to denote the type of mean size data.
#' \item MS_cv - The coefficient of variation of the observed mean size. If there are multiple fleets, a vector of length nfleet.
#' Default is 0.2.
#' \item s_CAA - Survey age composition data, an array of dimension nyears, maxage+1, nsurvey.
#' \item s_CAL - Survey length composition data, an array of dimension nyears, length(length_bin), nsurvey.
#' \item length_bin - A vector for the midpoints of the length bins for \code{CAL} and \code{s_CAL}. All bin widths should be equal in size.
#' \item C_eq - A numeric vector of length nfleet for the equilibrium catch for each fleet in \code{Chist} prior to the first year of the operating model.
#' Zero (default) implies unfished conditions in year one. Otherwise, this is used to estimate depletion in the first year of the data. Alternatively,
#' if one has a full CAA matrix, one could instead estimate "artificial" rec devs to generate the initial numbers-at-age (and hence initial depletion) in the first year of the model (see additional arguments).
#' \item C_eq_sd - A vector of standard deviations (lognormal distribution) for the equilibrium catches in \code{C_eq}.
#' If not provided, the default is 0.01. Only used if \code{condition = "catch"}.
#' \item E_eq - The equilibrium effort for each fleet in \code{Ehist} prior to the first year of the operating model.
#' Zero (default) implies unfished conditions in year one. Otherwise, this is used to estimate depletion in the first year of the data.
#' \item abs_I - Optional, an integer vector to indicate which indices are in absolute magnitude. Use 1 to set q = 1, otherwise use 0 to estimate q.
#' \item I_units - Optional, an integer vector to indicate whether indices are biomass based (1) or abundance-based (0). By default, all are biomass-based.
#' \item age_error - Optional, a square matrix of maxage + 1 rows and columns to specify ageing error. The aa-th column assigns a proportion of the true age in the
#' a-th row to observed age. Thus, all rows should sum to 1. Default is an identity matrix (no ageing error).
#' \item sel_block - Optional, for time-varying fleet selectivity (in time blocks), a integer matrix of nyears rows and nfleet columns to assigns a selectivity function to a fleet for certain years.
#' }
#'  
#' @section Additional arguments:
#' For \code{RCM}, additional arguments can be passed to the model via \code{...}:
#'
#' \itemize{
#' \item plusgroup: Logical for whether the maximum age is a plusgroup or not. By default, TRUE.
#' \item fix_dome: Logical for whether the dome selectivity parameter for fleets is fixed. Used primarily for backwards compatibility, this is overridden by the map argument.
#' \item resample: Logical, whether the OM conditioning parameters (recruitment, fishing mortality, SSB, selectivity, etc.) are obtained by sampling the Hessian matrix from
#' a single model fit. By default FALSE. This feature requires identical biological parameters among simulations.
#' }
#' @section start and map:
#' Starting values can be specified in a named list for the following:
#' \itemize{
#' \item vul_par: A matrix of 3 rows and nfleet columns for starting values for fleet selectivity. The three rows correspond
#' to LFS (length of full selectivity), L5 (length of 5 percent selectivity), and Vmaxlen (selectivity at length Linf). By default,
#' the starting values are values from the OM object. If any selectivity = "free", then this matrix needs to be of maxage+1 rows where
#' the row specifies the selectivity at age. See Articles section.
#' \item ivul_par: A matrix of 3 rows and nsurvey columns for starting values for fleet selectivity. Same setup as vul_par. Values in the column are ignored
#' if \code{s_selectivity} is mapped to a fishing fleet (add NA placeholders in that case). 
#' If any \code{s_selectivity = "free"}, then this matrix needs to be of maxage+1 rows where
#' the row specifies the selectivity at age. 
#' \item log_rec_dev: A numeric vector of length nyears for the starting values of the log-recruitment deviations.
#' \item log_early_rec_dev: A numeric vector of length OM@@maxage for the starting values of the recruitment deviations controlling the abundance-at-age in the first year of the model.
#' \item q: A numeric vector of length nsurvey for index catchability. See \href{https://openmse.com/tutorial-rcm-select/}{online article} for more information.
#' }
#' 
#' Parameters can be fixed with the map argument (also a named list, corresponding to the start list). Each
#' vector or matrix in the map argument will be the same dimension as in the start entry. If an entry is \code{NA}, the corresponding parameter is fixed in the model to the starting
#' value. Otherwise, an integer for each independent parameter, i.e., shared or mirrored parameters get the same integer entry.
#' 
#' \itemize{
#' \item vul_par: An integer matrix of the same dimension as vul_par. By default, selectivity is fixed if there are no age or length composition for that fleet
#' or survey, otherwise estimated. Unused cells in the vul_par matrix should be given NA in the map matrix.
#' \item ivul_par: The map argument for the survey selectivity parameters (same dimension as ivul_par). Placeholder parameters should have a map value of NA.
#' \item log_early_rec_dev: A vector of length OM@@maxage that indexes which recruitment deviates for the cohorts in the first year of the model are fixed (using NA) or estimated (a separate integer).
#' By default, no deviates are estimated (all are NA).
#' \item log_rec_dev: A vector of length OM@@nyears that indexes which recruitment deviates are fixed (using NA) or estimated (a separate integer).
#' By default, all these deviates are estimated.
#' \item q: A vector of length nsurvey for index catchability. q should be an estimated parameter when sharing across surveys (perhaps with differing selectivity). Otherwise, it is solved analytically
#' where individual parameters are independent of other indices. Use \code{RCMdata@abs_I} for fixing the catchability to 1. See \href{https://openmse.com/tutorial-rcm-select/}{online article} for more information.
#' }
#' 
#' @section Likelihood weights:
#' \code{LWT} is an optional named list containing the likelihood weights (values >= 0) with the possible options:
#' \itemize{
#' \item Chist, CAA, CAL, MS, C_eq: A vector of length nfleet for each.
#' \item Index, IAA, IAL: A vector of length nsurvey for each.
#' }
#'
#' By default, all likelihood weights are equal to one if not specified by the user.
#'
#' Annual multinomial sample sizes for the age and length comps can now be provided directly in the 
#' \linkS4class{RCMdata} object. For a list or \linkS4class{Data} object, use the \code{ESS} argument.
#' @author Q. Huynh
#' @examples 
#' \donttest{ 
#' # An example that conditions a Pacific cod operating model. There are 48 simulations, 
#' # where values of natural mortality and steepness are sampled from distributions. 
#' # The model is fitted with priors on the index catchability. Maturity and selectivity 
#' # are knife-edge at the age of 2 years. See online tutorial for more information.
#' 
#' data(pcod) 
#' mat_ogive <- pcod$OM@cpars$Mat_age[1, , 1]
#' out <- RCM(OM = pcod$OM, data = pcod$data, 
#'            condition = "catch", mean_fit = TRUE,
#'            selectivity = "free", s_selectivity = rep("SSB", ncol(pcod$data@Index)),
#'            start = list(vul_par = matrix(mat_ogive, length(mat_ogive), 1)),
#'            map = list(vul_par = matrix(NA, length(mat_ogive), 1),
#'                       log_early_rec_dev = rep(1, pcod$OM@maxage)),
#'            prior = pcod$prior)
#' plot(out, s_name = colnames(pcod$data@Index))
#' 
#' # Alternative OM with age-3 maturity and selectivity instead.
#' out_age3 <- local({
#'   pcod$OM@cpars$Mat_age[, 2, ] <- 0
#'   mat_ogive_age3 <- pcod$OM@cpars$Mat_age[1, , 1]
#'   RCM(OM = pcod$OM, data = pcod$data, 
#'       condition = "catch", mean_fit = TRUE,
#'       selectivity = "free", s_selectivity = rep("SSB", ncol(pcod$data@Index)),
#'       start = list(vul_par = matrix(mat_ogive_age3, length(mat_ogive_age3), 1)),
#'       map = list(vul_par = matrix(NA, length(mat_ogive_age3), 1),   
#'                  log_early_rec_dev = rep(1, pcod$OM@maxage)),
#'       prior = pcod$prior)
#' })
#'   
#' compare_RCM(out, out_age3, scenario = list(names = c("Age-2 maturity", "Age-3 maturity")),
#'             s_name = colnames(pcod$data@Index))
#'              
#' Hist <- runMSE(out@OM, Hist = TRUE)            
#' } 
#' 
#' @seealso \link{plot.RCModel} \linkS4class{RCModel} \link{compare_RCM} \link{pcod} \link{RCM2MOM} \link{posterior}
#' @references 
#' Thorson et al. 2017. Model-based estimates of effective sample size in stock assessment models using the Dirichlet-multinomial distribution.
#' Fish. Res. 192:84-93. \doi{10.1016/j.fishres.2016.06.005}
#' @importFrom dplyr %>%
#' @export
setGeneric("RCM", function(OM, data, ...) standardGeneric("RCM"))

#' @rdname RCM
#' @export
setMethod("RCM", signature(OM = "OM", data = "RCMdata"),
          function(OM, data, condition = "catch", selectivity = "logistic", s_selectivity = NULL, LWT = list(),
                   comp_like = c("multinomial", "lognormal", "mvlogistic", "dirmult1", "dirmult2"), prior = list(),
                   max_F = 3, cores = 1L, integrate = FALSE, mean_fit = FALSE, drop_nonconv = FALSE,
                   drop_highF = FALSE, control = list(iter.max = 2e+05, eval.max = 4e+05), 
                   start = list(), map = list(), silent = FALSE, ...) {
            RCM_int(OM = OM, RCMdata = data, condition = condition, selectivity = selectivity, s_selectivity = s_selectivity, LWT = LWT,
                    comp_like = comp_like, prior = prior, max_F = max_F, cores = cores, integrate = integrate, mean_fit = mean_fit,
                    drop_nonconv = drop_nonconv, drop_highF = drop_highF, control = control, start = start, map = map, silent = silent, ...)
          })

#' @rdname RCM
#' @export
setMethod("RCM", signature(OM = "OM", data = "list"),
          function(OM, data, condition = "catch", selectivity = "logistic", s_selectivity = NULL, LWT = list(),
                   comp_like = c("multinomial", "lognormal", "mvlogistic", "dirmult1", "dirmult2"), ESS = c(30, 30), prior = list(),
                   max_F = 3, cores = 1L, integrate = FALSE, mean_fit = FALSE, drop_nonconv = FALSE,
                   drop_highF = FALSE, control = list(iter.max = 2e+05, eval.max = 4e+05),
                   start = list(), map = list(), silent = FALSE, ...) {
            
            .Deprecated(msg = "Using a list of input data to RCM is now deprecated. Use an RCMdata object, i.e., new(\"RCMdata\")")
            
            dataS4 <- new("RCMdata")
            if (length(ESS) == 1) ESS <- rep(ESS, 2)
            
            if (!is.null(data$Chist)) dataS4@Chist <- data$Chist 
            if (!is.null(data$C_sd)) dataS4@C_sd <- data$C_sd 
            if (!is.null(data$Ehist)) dataS4@Ehist <- data$Ehist
            
            if (!is.null(data$CAA)) {
              dataS4@CAA <- data$CAA
              if (is.matrix(dataS4@CAA)) {
                dataS4@CAA_ESS <- apply(dataS4@CAA, 1, sum, na.rm = TRUE) %>% pmin(ESS[1])
              } else {
                dataS4@CAA_ESS <- apply(dataS4@CAA, c(1, 3), sum, na.rm = TRUE) %>% pmin(ESS[1])
              }
            }
            if (!is.null(data$CAL)) {
              dataS4@CAL <- data$CAL 
              if (is.matrix(dataS4@CAL)) {
                dataS4@CAL_ESS <- apply(dataS4@CAL, 1, sum, na.rm = TRUE) %>% pmin(ESS[2])
              } else {
                dataS4@CAL_ESS <- apply(dataS4@CAL, c(1, 3), sum, na.rm = TRUE) %>% pmin(ESS[2])
              }
            }
            
            if (!is.null(data$length_bin)) dataS4@length_bin <- data$length_bin 
            if (!is.null(data$MS)) dataS4@MS <- data$MS
            if (!is.null(data$MS_type)) dataS4@MS_type <- data$MS_type
            if (!is.null(data$MS_cv)) dataS4@MS_cv <- data$MS_cv
            
            if (!is.null(data$Index)) dataS4@Index <- data$Index
            if (!is.null(data$I_sd)) dataS4@I_sd <- data$I_sd
            
            if (!is.null(data$s_CAA)) dataS4@IAA <- data$s_CAA
            if (!is.null(data$IAA)) dataS4@IAA <- data$IAA
            if (length(dataS4@IAA)) {
              if (is.matrix(dataS4@IAA)) {
                dataS4@IAA_ESS <- apply(dataS4@IAA, 1, sum, na.rm = TRUE) %>% pmin(ESS[1])
              } else {
                dataS4@IAA_ESS <- apply(dataS4@IAA, c(1, 3), sum, na.rm = TRUE) %>% pmin(ESS[1])
              }
            }
            
            if (!is.null(data$s_CAL)) dataS4@IAL <- data$s_CAL
            if (!is.null(data$IAL)) dataS4@IAL <- data$IAL
            if (length(dataS4@IAL)) {
              if (is.matrix(dataS4@IAL)) {
                dataS4@IAL_ESS <- apply(dataS4@IAL, 1, sum, na.rm = TRUE) %>% pmin(ESS[2])
              } else {
                dataS4@IAL_ESS <- apply(dataS4@IAL, c(1, 3), sum, na.rm = TRUE) %>% pmin(ESS[2])
              }
            }
            
            if (!is.null(data$C_eq)) dataS4@C_eq <- data$C_eq
            if (!is.null(data$C_eq_sd)) dataS4@C_eq_sd <- data$C_eq_sd
            if (!is.null(data$E_eq)) dataS4@E_eq <- data$E_eq
            
            if (!is.null(data$abs_I)) dataS4@abs_I <- data$abs_I
            if (!is.null(data$I_units)) dataS4@I_units <- data$I_units
            
            if (!is.null(data$age_error)) dataS4@age_error <- data$age_error
            if (!is.null(data$sel_block)) dataS4@sel_block <- data$sel_block
            
            RCM_int(OM = OM, RCMdata = dataS4, condition = condition, selectivity = selectivity, s_selectivity = s_selectivity, LWT = LWT,
                    comp_like = comp_like, prior = prior, max_F = max_F, cores = cores, integrate = integrate, mean_fit = mean_fit,
                    drop_nonconv = drop_nonconv, drop_highF = drop_highF, control = control, start = start, map = map, silent = silent, ...)
          })


#' @rdname RCM
#' @export
setMethod("RCM", signature(OM = "OM", data = "Data"),
          function(OM, data, condition = "catch", selectivity = "logistic", s_selectivity = NULL, LWT = list(),
                   comp_like = c("multinomial", "lognormal", "mvlogistic", "dirmult1", "dirmult2"), ESS = c(30, 30), prior = list(),
                   max_F = 3, cores = 1L, integrate = FALSE, mean_fit = FALSE, drop_nonconv = FALSE,
                   drop_highF = FALSE, control = list(iter.max = 2e+05, eval.max = 4e+05), 
                   start = list(), map = list(), silent = FALSE, ...) {

            condition <- match.arg(condition, choices = c("catch", "catch2", "effort"))
            extra_args <- list(...)
            if (length(ESS) == 1) ESS <- rep(ESS, 2)
            
            ####### Check maxage from Data and OM
            if (!is.na(data@MaxAge) && OM@maxage != data@MaxAge) {
              warning("Data@MaxAge is not equal to OM@maxage")
            }

            ####### Catch and ML from Data object
            vec_slot <- c("Cat", "ML", "CV_Cat")
            data_vec <- lapply(vec_slot, vec_slot_fn, Data = data) %>% structure(names = vec_slot)

            ####### Age/length comps from Data object
            matrix_slot <- c("CAA", "CAL")
            data_matrix <- lapply(matrix_slot, matrix_slot_fn, Data = data) %>% structure(names = matrix_slot)

            ####### Generate data list for RCM
            dataS4 <- new("RCMdata")

            # Catch or effort
            if (condition == "effort") {
              if (all(is.na(data@Effort))) {
                stop("Conditioning on effort but no effort series found", call. = FALSE)
              } else {
                dataS4@Ehist <- data@Effort[1, ]
                if (nrow(data@Effort) == OM@nsim) { # Add sketched effort matrix to OM
                  OM@cpars$Find <- data@Effort
                  extra_args$OMeff <- TRUE
                }
              }
            } else if (is.null(data_vec$Cat)) {
              stop("Conditioning on catch but no catch data found", call. = FALSE)
            }
            if (!is.null(data_vec$Cat)) dataS4@Chist <- data_vec$Cat
            if (!is.null(data_vec$CV_Cat)) dataS4@C_sd <- sdconv(1, data_vec$CV_Cat) 

            # Index
            Ind <- pull_Index(data, OM@maxage)
            if (!is.null(Ind$Index)) {
              dataS4@Index <- Ind$Index
              dataS4@I_sd <- Ind$I_sd
              if (is.null(Ind$I_sd)) stop("No standard errors for the index was found.")

              if (any(!is.na(Ind$V))) {
                extra_args$ivul_par <- Ind$V
                extra_args$map_ivul_par <- array(NA, dim = dim(Ind$V))
              }
              dataS4@I_units <- Ind$I_units
            }

            # Length/age comps
            if (!is.null(data_matrix$CAA)) {
              dataS4@CAA <- data_matrix$CAA
              dataS4@CAA_ESS <- apply(dataS4@CAA, 1, sum, na.rm = TRUE) %>% pmin(ESS[1])
            }
            if (!is.null(data_matrix$CAL)) {
              dataS4@CAL <- data_matrix$CAL
              dataS4@CAL_ESS <- apply(dataS4@CAL, 1, sum, na.rm = TRUE) %>% pmin(ESS[2])
            }
            if (!is.null(data_matrix$CAL)) {
              if (all(is.na(data@CAL_mids))) {
                stop("No length bins found in Data@CAL_mids.", call. = FALSE)
              } else dataS4@length_bin <- data@CAL_mids
            }

            # Mean length
            if (!is.null(data_vec$ML)) dataS4@MS <- data_vec$ML # By default, MS_type = "length" and CV = 0.2

            # Equilibrium catches and/or effort - nothing happens if NULL
            if (!is.null(extra_args$C_eq)) dataS4@C_eq <- extra_args$C_eq
            if (!is.null(extra_args$E_eq)) dataS4@E_eq <- extra_args$E_eq
            
            if (!is.null(extra_args$ivul_par)) start$ivul_par <- extra_args$ivul_par
            if (!is.null(extra_args$map_ivul_par)) map$ivul_par <- extra_args$map_ivul_par

            ####### Run RCM
            output <- RCM_int(OM = OM, RCMdata = dataS4, condition = condition, selectivity = selectivity, s_selectivity = Ind$s_sel, LWT = LWT,
                              comp_like = comp_like, prior = prior, max_F = max_F, cores = cores, integrate = integrate, mean_fit = mean_fit,
                              drop_nonconv = drop_nonconv, drop_highF = drop_highF, control = control,
                              map = map, start = start, silent = silent, OMeff = extra_args$OMeff, ...)
            
            ####### Re-assign index slots from AddInd to their original places
            if (any(Ind$slotname != "AddInd")) {
              Data_out <- output@OM@cpars$Data
              for(i in 1:length(Ind$slotname)) {
                if (Ind$slotname[i] != "AddInd") {
                  slot(Data_out, Ind$slotname[i]) <- Data_out@AddInd[, i, ] %>% 
                    matrix(1, length(Data_out@Year))
                  slot(Data_out, paste0("CV_", Ind$slotname[i])) <- Data_out@CV_AddInd[, i, ] %>% 
                    matrix(1, length(Data_out@Year))
                  
                  if (Ind$slotname[i] == "Ind") output@OM@cpars$I_beta <- rep(1, output@OM@nsim)
                  if (Ind$slotname[i] == "VInd") output@OM@cpars$VI_beta <- rep(1, output@OM@nsim)
                  if (Ind$slotname[i] == "SpInd") output@OM@cpars$SpI_beta <- rep(1, output@OM@nsim)
                }
              }

              ind <- Ind$slotname == "AddInd"
              if (all(!ind)) {
                Data_out@AddInd <- Data_out@CV_AddInd <- Data_out@AddIndV <- array(NA, c(1, 1, 1))
                Data_out@AddIndType <- NA
                
                output@OM@cpars$AddIbeta <- NULL
              } else {
                Data_out@AddInd <- Data_out@AddInd[, ind, , drop = FALSE]
                Data_out@CV_AddInd <- Data_out@CV_AddInd[, ind, , drop = FALSE]
                Data_out@AddIndV <- Data_out@AddIndV[, ind, , drop = FALSE]
                Data_out@AddIndType <- Data_out@AddIndType[ind]
                
                output@OM@cpars$AddIbeta <- output@OM@cpars$AddIbeta[, ind, drop = FALSE]
              }
              output@OM@cpars$Data <- Data_out
            }

            ####### Done.
            return(output)
          })
