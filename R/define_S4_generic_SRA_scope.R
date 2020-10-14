

#' @name SRA_scope
#' @aliases SRA_scope SRA_scope,OM,list-method SRA_scope,OM,Data-method
#' @title Stock-reduction analysis (SRA) for conditioning operating models
#'
#' @description Intended for conditioning operating models for DLMtool. For data-limited stocks, this function can generate a range of potential depletion scenarios inferred from sparse data.
#' From a historical time series of total catch or effort, and potentially age/length compositions and multiple indices of abundance, the SRA returns a range of values for depletion, selectivity,
#' unfished recruitment (R0), historical fishing effort, and recruitment deviations for the operating model. This is done by sampling life history parameters
#' provided by the user and fitting a statistical catch-at-age model (with the predicted catch equal to the observed catch).
#' Alternatively one can do a single model fit and sample the covariance matrix to generate an operating model with uncertainty based on the model fit.
#' Either a full catch (conditioned on catch) or effort (conditioned on effort) time series is needed but missing data (as NAs) are allowed for all other data types.
#'
#' @param OM An object of class \linkS4class{OM} that specifies natural mortality (M), growth (Linf, K, t0, a, b), stock-recruitment relationship,
#' steepness, maturity parameters (L50 and L50_95), standard deviation of recruitment variability (Perr), as well as index uncertainty (Iobs).
#' @param data Data inputs formatted in a list object (preferred). Alternatively, \code{data} can be a \linkS4class{Data} S4 object. See Data section below.
#' @param condition String to indicate whether the SRA model is conditioned on "catch" (where F are estimated parameters), "catch2" (where F is solved internally using Newton's method),
#' or "effort".
#' @param selectivity A character vector of length nfleet to indicate \code{"logistic"}, \code{"dome"}, or \code{"free"} selectivity for each fleet in \code{Chist}.
#' If there is time-varying selectivity, this is a character vector of length nsel_block (see Data section below). "free" indicates independent selectivity parameters for each age,
#' and additional modifications for fixing selectivity parameters will likely be needed. See Additional arguments section.
#' @param s_selectivity A vector of length nsurvey to indicate the selectivity of the corresponding columns in \code{data$Index}. Use \code{"B"} for
#' total biomass, or \code{"SSB"} for spawning biomass (by default, "B" is used). Use numbers if the survey selectivity follows a fleet (corresponding to the columns in data$Chist, e.g., 1 = first fleet/column and so on).
#' If the survey selectivity is otherwise independent of anything else in the model, use \code{"logistic"}, \code{"dome"}, or \code{"free"} to specify the functional form of selectivity, and
#' see Additional arguments section for setup of survey selectivity parameters. See \href{../doc/SRA_scope_sel.html}{selectivity vignette} for more information.
#' @param LWT A named list of likelihood weights for the SRA model. See below.
#' @param comp_like A string indicating either \code{"multinomial"} (default) or \code{"lognormal"} distributions for the composition data.
#' @param ESS If \code{comp_like = "multinomial"}, a numeric vector of length two to cap the maximum effective samples size of the age and length compositions,
#' respectively, for the multinomial likelihood function. The effective sample size of an age or length composition sample is the minimum of ESS or the number of observations
#' (sum across columns). For more flexibility, set ESS to be very large and alter the age and length arrays as needed.
#' @param max_F The maximum F for any fleet in the scoping model (higher F's in the model are penalized in the objective function). See also `drop_highF`.
#' @param cores Integer for the number of CPU cores for the stock reduction analysis.
#' @param integrate Logical, whether to treat recruitment deviations as penalized parameters in the likelihood (FALSE) or random effects to be marginalized out of the likelihood (TRUE).
#' @param mean_fit Logical, whether to run an additional with mean values of life history parameters from the OM.
#' @param sims A logical vector of length \code{OM@@nsim} or a numeric vector indicating which simulations to keep.
#' @param drop_nonconv Logical, whether to drop non-converged fits of the SRA model, including fits where F = NA.
#' @param drop_highF Logical, whether to drop fits of the SRA model where F = \code{max_F}.
#' @param control A named list of arguments (e.g, max. iterations, etc.) for optimization, to be passed to the control argument of \code{\link[stats]{nlminb}}.
#' @param ... Other arguments to pass in for starting values of parameters and fixing parameters. See details.
#'
#' @details
#' Fleet selectivity is fixed to values sampled from \code{OM} if no age or length compositions are provided.
#'
#' Survey selectivity is estimable only if \code{s_CAA} or \code{s_CAL} is provided. Otherwise, the selectivity should
#' be mirrored to a fleet (vulnerable biomass selectivity) or indexed to total or spawning biomass (see \code{s_selectivity}).
#'
#' Parameters that were used in the fitting model are placed in the \code{SRA@@OM@@cpars} list.

#' If the operating model \code{OM} uses time-varying growth or M, then those trends will be used in the SRA as well.
#' Time-varying life history parameters can create ambiguity in the calculation and interpretation of depletion and reference points in \link[DLMtool]{runMSE}.
#' See section D.5 of \code{DLMtool::userguide()}.
#'
#' The easiest way to turn off time-varying growth/M is by setting: \code{OM@@Msd <- OM@@Linfsd <- OM@@Ksd <- c(0, 0)}.
#'
#' \code{Sub_cpars} is a convenient function to subset simulations
#' for the operating model, for example, to remove simulations from unconverged model fits or outlier simulations.
#'
#' To play with alternative fits by excluding indices, for example, or other optional data, set the corresponding likelihood weight to zero. The model will still generate the inferred
#' index but the data won't enter the likelihood. See section on likelihood weights.
#'
#' @return An object of class \linkS4class{SRA} (see link for description of output).
#'
#' @section Vignette:
#' Three vignettes are available for the SRA model:
#'
#' \itemize{
#' \item \href{../doc/SRA_scope.html}{General overview of approach}
#' \item \href{../doc/SRA_scope_eq.html}{Mathematical description}
#' \item \href{../doc/SRA_scope_sel.html}{Setup of selectivity settings} (useful for more data-rich cases)
#' }
#'
#' @section Data:
#' One of indices, age compositions, or length compositions should be provided in addition to the historical catch or effort. Not all arguments
#' are needed to run the model (some have defaults, while others are ignored if not applicable depending on the data provided).
#'
#' The \code{data} variable can be a named list that includes:
#'
#' \itemize{
#' \item Chist - A vector of historical catch, should be of length OM@@nyears. If there are multiple fleets: a matrix of OM@@nyears rows and nfleet columns.
#' Ideally, the first year of the catch series represents unfished conditions (see also \code{C_eq}).
#' \item Ehist - A vector of historical effort, should be of length OM@@nyears (see also \code{E_eq}).
#' \item Index - A vector of values of an index (of length OM@@nyears). If there are multiple surveys: a matrix of historical indices of abundances, with rows
#' indexing years and columns indexing surveys. Age-specific indices should be numbers-specific while all others are weight-based.
#' \item I_sd - A vector or matrix of standard deviations (lognormal distribution) for the indices corresponding to the entries in \code{Index}.
#' If not provided, this function will use values from \code{OM@@Iobs}.
#' \item I_type - Obsolete as of version 2.0. See \code{s_selectivity} argument.
#' \item CAA - Fishery age composition matrix with nyears rows and OM@@maxage columns. If multiple fleets: an array with dimension: nyears, OM@@maxage, and nfleets.
#' \item CAL - Fishery length composition matrix with nyears rows and columns indexing the length bin. If multiple fleets: an array with dimension: nyears,
#' length bins, and nfleets.
#' \item MS - A vector of fishery mean size (MS, either mean length or mean weight) observations (length OM@@nyears), or if multiple fleets: matrix of dimension: nyears and nfleets.
#' Generally, mean lengths should not be used if \code{CAL} is also provided, unless mean length and length comps are independently sampled.
#' \item MS_type - A character (either \code{"length"} (default) or \code{"weight"}) to denote the type of mean size data.
#' \item MS_cv - The coefficient of variation of the observed mean size. If there are multiple fleets, a vector of length nfleet.
#' Default is 0.2.
#' \item s_CAA - Survey age composition data, an array of dimension nyears, maxage, nsurvey.
#' \item s_CAL - Survey length composition data, an array of dimension nyears, length(length_bin), nsurvey.
#' \item length_bin - A vector for the midpoints of the length bins for \code{CAL} and \code{s_CAL}. All bin widths should be equal in size.
#' \item C_eq - A numeric vector of length nfleet for the equilibrium catch for each fleet in \code{Chist} prior to the first year of the operating model.
#' Zero (default) implies unfished conditions in year one. Otherwise, this is used to estimate depletion in the first year of the data. Alternatively,
#' if one has a full CAA matrix, one could instead estimate "artificial" rec devs to generate the initial numbers-at-age (and hence initial depletion) in the first year of the model (see additional arguments).
#' \item E_eq - The equilibrium effort for each fleet in \code{Ehist} prior to the first year of the operating model.
#' Zero (default) implies unfished conditions in year one. Otherwise, this is used to estimate depletion in the first year of the data.
#' \item abs_I - Optional, an integer vector to indicate which indices are in absolute magnitude. Use 1 to set q = 1, otherwise use 0 to estimate q.
#' \item I_units - Optional, an integer vector to indicate whether indices are biomass based (1) or abundance-based (0). By default, all are biomass-based.
#' \item age_error - Optional, a square matrix of maxage rows and columns to specify ageing error. The aa-th column assigns a proportion of the true age in the
#' a-th row to observed age. Thus, all rows should sum to 1. Default is an identity matrix (no ageing error).
#' \item sel_block - Optional, for time-varying fleet selectivity (in time blocks), a integer matrix of nyears rows and nfleet columns to assigns a selectivity function to a fleet for certain years.
#' See the \href{../doc/SRA_scope_sel.html}{selectivity} vignette for more details.
#' }
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
#' There is no slot in the Data S4 object for the equilibrium catch/effort. These can be passed in the function call, i.e., \code{SRA_scope(OM, Data, C_eq = C_eq, ...)}.
#'
#'
#' @section Additional arguments:
#' For \code{SRA_scope}, additional arguments can be passed to the model via \code{...}:
#'
#' \itemize{
#' \item vul_par: A matrix of 3 rows and nfleet columns for starting values for fleet selectivity. The three rows correspond
#' to LFS (length of full selectivity), L5 (length of 5 percent selectivity), and Vmaxlen (selectivity at length Linf). By default,
#' the starting values are values from the OM object. If any selectivity = "free", then this matrix needs to be of maxage rows where
#' the row specifies the selectivity at age. See the \href{../doc/SRA_scope_sel.html}{selectivity} vignette for more information.
#' \item s_vul_par: A matrix of 3 rows and nsurvey columns for starting values for fleet selectivity. Same setup as vul_par. These values are only
#' used if \code{s_selectivity = "est"} for the corresponding fleet. Otherwise, placeholders should be used to complete the matrix.
#' \item map_vul_par: An integer matrix of the same dimension as vul_par. This is the 'map' argument for vul_par in TMB, see \link[TMB]{MakeADFun}, which indicates whether selectivity parameters are fixed
#' or estimated. If an entry is \code{NA}, the corresponding parameter is fixed in the model to the starting
#' value. Otherwise, an integer for each independent parameter. By default, selectivity is fixed if there are no age or length composition for that fleet
#' or survey, otherwise estimated. Unused cells in the vul_par matrix should be given NA in the map matrix.
#' \item map_s_vul_par: The map argument for the survey selectivity parameters (same dimension as s_vul_par). Placeholder parameters should have a map value of NA.
#' \item map_log_early_rec_dev: A vector of length OM@@maxage - 1 that indexes which recruitment deviates for the cohorts in the first year of the model are fixed (using NA) or estimated (a separate integer).
#' By default, no deviates are estimated.
#' \item map_log_rec_dev: A vector of length OM@@nyears that indexes which recruitment deviates are fixed (using NA) or estimated (a separate integer).
#' By default, all deviates are estimated.
#' \item plusgroup: Logical for whether the maximum age is a plusgroup or not. By default, TRUE.
#' \item fix_dome: Logical for whether the dome selectivity parameter for fleets is fixed. Used primarily for backwards compatibility, this is overridden by map_vul_par.
#' \item resample: Logical, whether the OM conditioning parameters (recruitment, fishing mortality, SSB, selectivity, etc.) are obtained by sampling the Hessian matrix from
#' a single model fit. By default FALSE. This feature requires identical biological parameters among simulations.
#' }
#'
#' @section Likelihood weights:
#' \code{LWT} is an optional named list containing the likelihood weights (values >= 0) with the possible options:
#' \itemize{
#' \item Chist, CAA, CAL, MS, C_eq: A vector of length nfleet for each.
#' \item Index, s_CAA, s_CAL: A vector of length nsurvey for each.
#' }
#'
#' By default, all likelihood weights are equal to one if not specified by the user.
#'
#' Weighting for CAA and CAL can also be adjusted by changing the multinomial sample size. For \code{CAA}, \code{CAL}, \code{s_CAA}, and \code{s_CAL}, the arrays should be set up so that
#' the annual number of observations will be equal to the presumed multinomial sample size. Argument \code{ESS} provides a shortcut
#' to cap the multinomial sample size for age and length comps.
#'
#' @author Q. Huynh
#' @seealso \link{plot.SRA} \linkS4class{SRA}
#' @importFrom dplyr %>%
#' @export
setGeneric("SRA_scope", function(OM, data, ...) standardGeneric("SRA_scope"))

#' @rdname SRA_scope
#' @export
setMethod("SRA_scope", signature(OM = "OM", data = "list"),
          function(OM, data, condition = c("catch", "catch2", "effort"), selectivity = "logistic", s_selectivity = NULL, LWT = list(),
                   comp_like = c("multinomial", "lognormal"), ESS = c(30, 30),
                   max_F = 3, cores = 1L, integrate = FALSE, mean_fit = FALSE, drop_nonconv = FALSE,
                   drop_highF = FALSE, control = list(iter.max = 2e+05, eval.max = 4e+05), ...) {

            SRA_scope_int(OM = OM, data = data, condition = condition, selectivity = selectivity, s_selectivity = s_selectivity, LWT = LWT,
                          comp_like = comp_like, ESS = ESS, max_F = max_F, cores = cores, integrate = integrate, mean_fit = mean_fit,
                          drop_nonconv = drop_nonconv, drop_highF = drop_highF, control = control, ...)

          })


#' @rdname SRA_scope
#' @export
setMethod("SRA_scope", signature(OM = "OM", data = "Data"),
          function(OM, data, condition = c("catch", "catch2", "effort"), selectivity = "logistic", s_selectivity = NULL, LWT = list(),
                   comp_like = c("multinomial", "lognormal"), ESS = c(30, 30),
                   max_F = 3, cores = 1L, integrate = FALSE, mean_fit = FALSE, drop_nonconv = FALSE,
                   drop_highF = FALSE, control = list(iter.max = 2e+05, eval.max = 4e+05), ...) {

            condition <- match.arg(condition)
            extra_args <- list(...)

            ####### Harvest catch and ML from Data object
            vec_slot <- c("Cat", "ML")
            data_vec <- lapply(vec_slot, vec_slot_fn, Data = data) %>% structure(names = vec_slot)

            ####### Harvest age/length comps from Data object
            matrix_slot <- c("CAA", "CAL")
            data_matrix <- lapply(matrix_slot, matrix_slot_fn, Data = data) %>% structure(names = matrix_slot)

            ####### Generate data list for SRA_scope
            data_list <- list()

            # Catch or effort
            if(condition == "effort") {
              if(all(is.na(data@Effort))) {
                stop("Conditioning on effort but no effort series found", call. = FALSE)
              } else {
                data_list$Ehist <- data@Effort[1, ]
                if(nrow(data@Effort) == OM@nsim) { # Add sketched effort matrix to OM
                  OM@cpars$Find <- data@Effort
                  extra_args$OMeff <- TRUE
                }
              }
            } else {
              if(is.null(data_vec$Cat)) {
                stop("Conditioning on catch but no catch data found", call. = FALSE)
              } else {
                data_list$Chist <- data_vec$Cat
              }
            }

            # Index
            Ind <- pull_Index(data, OM@maxage)
            if(!is.null(Ind$Index)) {
              data_list$Index <- Ind$Index
              data_list$I_sd <- Ind$I_sd

              if(any(!is.na(Ind$V))) {
                extra_args$s_vul_par <- Ind$V
                extra_args$map_s_vul_par <- array(NA, dim = dim(Ind$V))
              }
              data_list$I_units <- Ind$I_units
            }

            # Length/age comps
            data_list$CAA <- data_matrix$CAA
            data_list$CAL <- data_matrix$CAL
            if(!is.null(data_list$CAL)) {
              if(all(is.na(data@CAL_mids))) {
                stop("No length bins found in Data@CAL_mids.", call. = FALSE)
              } else data_list$length_bin <- data@CAL_mids
            }

            # Mean length
            data_list$MS <- data_vec$ML # By default, MS_type = "length" and CV = 0.2

            # Equilibrium catches and/or effort - nothing happens if NULL
            data_list$C_eq <- extra_args$C_eq
            data_list$E_eq <- extra_args$E_eq

            ####### Run SRA_scope
            output <-
              SRA_scope_int(OM = OM, data = data_list, condition = condition, selectivity = selectivity, s_selectivity = Ind$s_sel, LWT = LWT,
                            comp_like = comp_like, ESS = ESS, max_F = max_F, cores = cores, integrate = integrate, mean_fit = mean_fit,
                            drop_nonconv = drop_nonconv, drop_highF = drop_highF, control = control,
                            OMeff = extra_args$OMeff, s_vul_par = extra_args$s_vul_par, map_s_vul_par = extra_args$map_s_vul_par, ...)

            ####### Re-assign index slots from AddInd to their original places
            if(any(Ind$slotname != "AddInd")) {
              Data_out <- output@OM@cpars$Data
              for(i in 1:length(Ind$slotname)) {
                if(Ind$slotname[i] != "AddInd") {
                  slot(Data_out, Ind$slotname[i]) <- Data_out@AddInd[, i, ]
                  slot(Data_out, paste0("CV_", Ind$slotname[i])) <- Data_out@CV_AddInd[, i, ]
                }
              }

              ind <- Ind$slotname == "AddInd"
              if(all(!ind)) {
                Data_out@AddInd <- Data_out@CV_AddInd <- Data_out@AddIndV <- array(NA, c(1, 1, 1))
                Data_out@AddIndType <- NA
              } else {
                Data_out@AddInd <- Data_out@AddInd[, ind, , drop = FALSE]
                Data_out@CV_AddInd <- Data_out@CV_AddInd[, ind, , drop = FALSE]
                Data_out@AddIndV <- Data_out@AddIndV[, ind, , drop = FALSE]
                Data_out@AddIndType <- Data_out@AddIndType[ind]
              }
              output@OM@cpars$Data <- Data_out
            }

            ####### Done.
            return(output)
          })
