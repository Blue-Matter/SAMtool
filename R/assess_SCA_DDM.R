
#' SCA model with density-dependent M
#' 
#' A modification of \code{SCA} that incorporates density-dependent effects on M based on biomass depletion (Forrest et al. 2018).
#' Set the bounds of M in the \code{M_bounds} argument, a length-2 vector where the first entry is M0, the M as B/B0 >= 1,
#' and the second entry is M1, the M as B/B0 approaches zero. Note that M0 can be greater than M1 (compensatory) 
#' or M0 can be less than M1 (depensatory).
#' 
#' @param x A position in the Data object (by default, equal to one for assessments).
#' @param Data An object of class Data
#' @param AddInd A vector of integers or character strings indicating the indices to be used in the model. Integers assign the index to
#' the corresponding index in Data@@AddInd, "B" (or 0) represents total biomass in Data@@Ind, "VB" represents vulnerable biomass in
#' Data@@VInd, and "SSB" represents spawning stock biomass in Data@@SpInd. Vulnerability to the survey is fixed in the model.
#' @param SR Stock-recruit function (either \code{"BH"} for Beverton-Holt, \code{"Ricker"}, or \code{"none"} for constant mean recruitment).
#' @param vulnerability Whether estimated vulnerability is \code{"logistic"} or \code{"dome"} (double-normal).
#' See details for parameterization.
#' @param catch_eq Whether to use the Baranov equation or Pope's approximation to calculate the predicted catch at age in the model.
#' @param CAA_dist Whether a multinomial or lognormal distribution is used for likelihood of the catch-at-age matrix. See details.
#' @param CAA_multiplier Numeric for data weighting of catch-at-age matrix if \code{CAA_hist = "multinomial"}. Otherwise ignored. See details.
#' @param rescale A multiplicative factor that rescales the catch in the assessment model, which
#' can improve convergence. By default, \code{"mean1"} scales the catch so that time series mean is 1, otherwise a numeric.
#' Output is re-converted back to original units.
#' @param max_age Integer, the maximum age (plus-group) in the model.
#' @param start Optional list of starting values. Entries can be expressions that are evaluated in the function. See details.
#' @param prior A named list (R0, h, M, and q) to provide the mean and standard deviations of prior distributions for those parameters. R0, index q, and M priors are
#' lognormal (provide the mean in normal space, SD in lognormal space). Beverton-Holt steepness uses a beta prior, while Ricker steepness uses a normal prior.
#' For index q, provide a matrix for nsurvey rows and 2 columns (for mean and SD), with NA in rows corresponding to indices without priors
#' For all others, provide a length-2 vector for the mean and SD.
#' See vignette for full description.
#' @param fix_h Logical, whether to fix steepness to value in \code{Data@@steep} in the model for \code{SCA}. This only affects
#' calculation of reference points for \code{SCA2}.
#' @param fix_F_equilibrium Logical, whether the equilibrium fishing mortality prior to the first year of the model
#' is estimated. If \code{TRUE}, \code{F_equilibrium} is fixed to value provided in \code{start} (if provided),
#' otherwise, equal to zero (assumes unfished conditions).
#' @param fix_omega Logical, whether the standard deviation of the catch is fixed. If \code{TRUE},
#' omega is fixed to value provided in \code{start} (if provided), otherwise, value based on \code{Data@@CV_Cat}.
#' @param fix_tau Logical, the standard deviation of the recruitment deviations is fixed. If \code{TRUE},
#' tau is fixed to value provided in \code{start} (if provided), otherwise, value based on \code{Data@@sigmaR}.
#' @param LWT A vector of likelihood weights for each survey.
#' @param common_dev Typically, a numeric for the number of most recent years in which a common recruitment deviation will
#' be estimated (in \code{SCA2}, uninformative years will have a recruitment closer to the mean, which can be very misleading,
#' especially near the end of the time series). By default, \code{"comp50"} uses the number of ages (smaller than the mode)
#' for which the catch-at-age matrix has less than half the abundance than that at the mode.
#' @param early_dev Numeric or character string describing the years for which recruitment deviations are estimated in \code{SCA}. By default,
#' equal to \code{"comp_onegen"}, where rec devs are estimated one full generation prior to the first year when catch-at-age (CAA) data are available.
#' With \code{"comp"}, rec devs are estimated starting in the first year with CAA. With \code{"all"}, rec devs start at the beginning of the model.
#' If numeric, the number of years after the first year of the model for which to start estimating rec devs. Use negative numbers for years prior to the first year.
#' @param late_dev Typically, a numeric for the number of most recent years in which recruitment deviations will
#' not be estimated in \code{SCA} (recruitment in these years will be based on the mean predicted by stock-recruit relationship).
#' By default, \code{"comp50"} uses the number of ages (smaller than the mode)
#' for which the catch-at-age matrix has less than half the abundance than that at the mode.
#' @param M_bounds A numeric vector of length 2 to indicate the M as B/B0 approaches zero and one, respectively.
#' By default, set to 75\% and 125\%, respectively, of Data@@Mort[x].
#' @param integrate Logical, whether the likelihood of the model integrates over the likelihood
#' of the recruitment deviations (thus, treating it as a random effects/state-space variable).
#' Otherwise, recruitment deviations are penalized parameters.
#' @param silent Logical, passed to \code{\link[TMB]{MakeADFun}}, whether TMB
#' will print trace information during optimization. Used for diagnostics for model convergence.
#' @param opt_hess Logical, whether the hessian function will be passed to \code{\link[stats]{nlminb}} during optimization
#' (this generally reduces the number of iterations to convergence, but is memory and time intensive and does not guarantee an increase
#' in convergence rate). Ignored if \code{integrate = TRUE}.
#' @param n_restart The number of restarts (calls to \code{\link[stats]{nlminb}}) in the optimization procedure, so long as the model
#' hasn't converged. The optimization continues from the parameters from the previous (re)start.
#' @param control A named list of arguments for optimization to be passed to
#' \code{\link[stats]{nlminb}}.
#' @param inner.control A named list of arguments for optimization of the random effects, which
#' is passed on to \code{\link[TMB]{newton}}.
#' @param ... Other arguments to be passed.
#' @details
#' See \link{SCA} for more information on all arguments.
#' @references
#' Forrest, R.E., Holt, K.R., and Kronlund, A.R. 2018. Performance of alternative harvest control rules for two Pacific groundfish
#' stocks with uncertain natural mortality: Bias, robustness and trade-offs. Fisheries Research 2016: 259-286.
#'
#' @examples
#' res <- SCA_DDM(Data = MSEtool::SimulatedData)
#' 
#' @author Q. Huynh
#' @return An object of class \linkS4class{Assessment}.
#' @seealso \link{plot.Assessment} \link{summary.Assessment} \link{retrospective} \link{profile} \link{make_MP}
#' @export
SCA_DDM <- function(x = 1, Data, AddInd = "B", SR = c("BH", "Ricker", "none"), 
                    vulnerability = c("logistic", "dome"), catch_eq = c("Baranov", "Pope"),
                    CAA_dist = c("multinomial", "lognormal"),
                    CAA_multiplier = 50, rescale = "mean1", max_age = Data@MaxAge,
                    start = NULL, prior = list(), fix_h = TRUE, fix_F_equilibrium = TRUE, fix_omega = TRUE, fix_tau = TRUE,
                    LWT = NULL, early_dev = c("comp_onegen", "comp", "all"), late_dev = "comp50", M_bounds = NULL,
                    integrate = FALSE, silent = TRUE, opt_hess = FALSE, n_restart = ifelse(opt_hess, 0, 1),
                    control = list(iter.max = 2e5, eval.max = 4e5), inner.control = list(), ...) {
  
  out <- SCA_(x, Data, AddInd, SR, vulnerability, catch_eq, CAA_dist, CAA_multiplier, rescale, max_age,
              start, prior, fix_h, fix_F_equilibrium, fix_omega, fix_tau, 
              tv_M = "DD", M_bounds = M_bounds,
              LWT = LWT, early_dev = early_dev, late_dev = late_dev, integrate = integrate,
              silent = silent, opt_hess = opt_hess, n_restart = n_restart,
              control = control, inner.control = inner.control, ...) 
  return(out)
}
class(SCA_DDM) <- "Assess"

