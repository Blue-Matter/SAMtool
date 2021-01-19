

#' @import MSEtool
#' @import TMB
#' @import graphics
#' @import stats
#' @import utils
#' @import methods
setOldClass("sdreport")
setClassUnion("sdreportAssess", members = c("character", "sdreport"))
setClassUnion("optAssess", members = c("list", "character"))
setClassUnion("vectormatrix", members = c("vector", "matrix"))



#' Class-\code{Assessment}
#'
#' An S4 class that contains assessment output. Created from a function of class \code{Assess}.
#'
#' @name Assessment-class
#' @docType class
#'
#' @slot Model Name of the assessment model.
#' @slot Name Name of Data object.
#' @slot conv Logical. Whether the assessment model converged (defined by whether TMB returned a
#' positive-definite covariance matrix for the model).
#' @slot UMSY Estimate of exploitation at maximum sustainable yield.
#' @slot FMSY Estimate of instantaneous fishing mortality rate at maximum sustainable yield.
#' @slot MSY Estimate of maximum sustainable yield.
#' @slot BMSY Biomass at maximum sustainable yield.
#' @slot SSBMSY Spawning stock biomass at maximum sustainable yield.
#' @slot VBMSY Vulnerable biomass at maximum sustainable yield.
#' @slot B0 Biomass at unfished equilibrium.
#' @slot R0 Recruitment at unfished equilibrium.
#' @slot N0 Abundance at unfished equilibrium.
#' @slot SSB0 Spawning stock biomass at unfished equilibrium.
#' @slot VB0 Vulnerable biomass at unfished equilibrium.
#' @slot h Steepness.
#' @slot U Time series of exploitation.
#' @slot U_UMSY Time series of relative exploitation.
#' @slot FMort Time series of instantaneous fishing mortality.
#' @slot F_FMSY Time series of fishing mortality relative to MSY.
#' @slot B Time series of biomass.
#' @slot B_BMSY Time series of biomass relative to MSY.
#' @slot B_B0 Time series of depletion.
#' @slot SSB Time series of spawning stock biomass.
#' @slot SSB_SSBMSY Time series of spawning stock biomass relative to MSY.
#' @slot SSB_SSB0 Time series of spawning stock depletion.
#' @slot VB Time series of vulnerable biomass.
#' @slot VB_VBMSY Time series of vulnerable biomass relative to MSY.
#' @slot VB_VB0 Time series of vulnerable biomass depletion.
#' @slot R Time series of recruitment.
#' @slot N Time series of population abundance.
#' @slot N_at_age Time series of numbers-at-age matrix.
#' @slot Selectivity Selectivity-at-age matrix.
#' @slot Obs_Catch Observed catch.
#' @slot Obs_Index Observed index.
#' @slot Obs_C_at_age Observed catch-at-age matrix.
#' @slot Catch Predicted catch.
#' @slot Index Predicted index.
#' @slot C_at_age Predicted catch-at-age matrix.
#' @slot Dev A vector of estimated deviation parameters.
#' @slot Dev_type A description of the deviation parameters, e.g. "log recruitment deviations".
#' @slot NLL Negative log-likelihood. A vector for the total likelihood, integrated across random effects if applicable, components,
#' and penalty term (applied when \code{U > 0.975} in any year).
#' @slot SE_UMSY Standard error of UMSY estimate.
#' @slot SE_FMSY Standard error of FMSY estimate.
#' @slot SE_MSY Standard error of MSY estimate.
#' @slot SE_U_UMSY Standard error of U/UMSY.
#' @slot SE_F_FMSY Standard error of F/FMSY.
#' @slot SE_B_BMSY Standard error of B/BMSY.
#' @slot SE_B_B0 Standard error of B/B0.
#' @slot SE_SSB_SSBMSY Standard error of SSB/SSBMSY.
#' @slot SE_SSB_SSB0 Standard error of SSB/SSB0.
#' @slot SE_VB_VBMSY Standard error of VB/VBMSY.
#' @slot SE_VB_VB0 Standard error of VB/VB0.
#' @slot SE_Dev A vector of standard errors of the deviation parameters.
#' @slot info A list containing the data and starting values of estimated parameters
#' for the assessment.
#' @slot forecast A list containing components for forecasting:
#' \itemize{
#' \item \code{per_recruit} A data frame of SPR (spawning potential ratio) and YPR (yield-per-recruit), calculated for
#' a range of harvest rate of 0 - 0.99 or instantaneous F from 0 - 2.5 FMSY. 
#' \item \code{catch_eq} A function that calculates the catch for the next year (after the model terminal year) when an
#' apical F is provided. 
#' }
#' @slot obj A list with components returned from \code{\link[TMB]{MakeADFun}}.
#' @slot opt A list with components from calling \code{\link[stats]{nlminb}} to \code{obj}.
#' @slot SD A list (class sdreport) with parameter estimates and their standard errors, obtained from
#' \code{\link[TMB]{sdreport}}.
#' @slot TMB_report A list of model output reported from the TMB executable, i.e. \code{obj$report()}, and derived quantities (e.g. MSY).
#' @slot dependencies A character string of data types required for the assessment.
#' @examples
#' \donttest{
#' output <- DD_TMB(Data = MSEtool::SimulatedData)
#' class(output)
#' }
#' @seealso \link{plot.Assessment} \link{summary.Assessment} \link{retrospective} \link{profile} \link{make_MP}
#' @author Q. Huynh
#' @export Assessment
#' @exportClass Assessment
Assessment <- setClass("Assessment",
                       slots = c(Model = "character", Name = "character", conv = "logical", UMSY = "numeric", FMSY = "numeric",
                                 MSY = "numeric", BMSY = "numeric", SSBMSY = "numeric", VBMSY = "numeric",
                                 B0 = "numeric", R0 = "numeric", N0 = "numeric", SSB0 = "numeric", VB0 = "numeric",
                                 h = "numeric", U = "numeric", U_UMSY = "numeric", FMort = "numeric", F_FMSY  = "numeric",
                                 B = "numeric", B_BMSY = "numeric", B_B0 = "numeric",
                                 SSB = "numeric", SSB_SSBMSY = "numeric", SSB_SSB0 = "numeric", VB = "numeric",
                                 VB_VBMSY = "numeric", VB_VB0 = "numeric",
                                 R = "numeric", N = "numeric", N_at_age = "matrix",
                                 Selectivity = "matrix", Obs_Catch = "numeric", Obs_Index = "vectormatrix",
                                 Obs_C_at_age = "matrix", Catch = "numeric", Index = "vectormatrix",
                                 C_at_age = "matrix", Dev = "numeric", Dev_type = "character",
                                 NLL = "numeric", SE_UMSY = "numeric", SE_FMSY = "numeric", SE_MSY = "numeric",
                                 SE_U_UMSY = "numeric", SE_F_FMSY = "numeric",
                                 SE_B_BMSY = "numeric", SE_B_B0 = "numeric",
                                 SE_SSB_SSBMSY = "numeric", SE_SSB_SSB0 = "numeric",
                                 SE_VB_VBMSY = "numeric", SE_VB_VB0 = "numeric",
                                 SE_Dev = "numeric", info = "ANY", forecast = "list",
                                 obj = "list", opt = "optAssess", SD = "sdreportAssess",
                                 TMB_report = "list", dependencies = "character"))

#' @name summary.Assessment
#' @title Summary of Assessment object
#' @aliases summary,Assessment-method
#' @description Returns a summary of parameter estimates and output from an \linkS4class{Assessment} object.
#' @param object An object of class \linkS4class{Assessment}
#' @return A list of parameters.
#' @examples
#' output <- DD_TMB(Data = MSEtool::SimulatedData)
#' summary(output)
#' @exportMethod summary
setMethod("summary", signature(object = "Assessment"), function(object) {
  if(is.character(object@opt) || is.character(object@SD)) warning("Did model converge?")
  f <- get(paste0("summary_", object@Model))
  f(object)
})


#setMethod("show", signature(object = "Assessment"), function(object) print(summary(object)))

#' @name plot.Assessment
#' @aliases plot,Assessment,missing-method
#' @title Plot Assessment object
#' @description Produces HTML file (via markdown) figures of parameter estimates and output from an \linkS4class{Assessment} object.
#'
#' @param x An object of class \linkS4class{Assessment}.
#' @param y An object of class \linkS4class{retro}.
#' @param filename Character string for the name of the markdown and HTML files.
#' @param dir The directory in which the markdown and HTML files will be saved.
#' @param ret_yr If greater than zero, then a retrospective analysis will be performed and results will be reported. The integer here corresponds
#' to the number of peels (the maximum number of terminal years for which the data are removed).
#' @param open_file Logical, whether the HTML document is opened after it is rendered.
#' @param quiet Logical, whether to silence the markdown rendering function.
#' @param render_args Arguments to pass to \link[rmarkdown]{render}.
#' @param ... Other arguments.
#' @return Returns invisibly the output from \link[rmarkdown]{render}.
#' @examples
#' output <- DD_TMB(Data = Simulation_1)
#' 
#' \donttest{
#' plot(output)
#' }
#' @seealso \link{retrospective}
#' @exportMethod plot
setMethod("plot", signature(x = "Assessment", y = "missing"),
          function(x, filename = paste0("report_", x@Model), dir = tempdir(), ret_yr = 0L,
                   open_file = TRUE, quiet = TRUE, render_args = list(), ...) {
            # Generating retrospective
            if(is.numeric(ret_yr) && ret_yr > 0) {
              ret_yr <- as.integer(ret_yr)
              message("Running retrospective...")
              ret <- retrospective(x, nyr = ret_yr, figure = FALSE)
            } else ret <- NULL
            report(x, ret, filename = filename, dir = dir, open_file = open_file, quiet = quiet, render_args = render_args, ...)
          })


if(getRversion() >= "2.15.1") {
  # Define global variables for Assessment objects
  utils::globalVariables(slotNames("Assessment"))
}


