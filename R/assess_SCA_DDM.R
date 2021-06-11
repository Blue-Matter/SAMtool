
#' SCA models with time-varying natural mortality
#' 
#' A modification of \code{SCA} that incorporates density-dependent effects on M based on biomass depletion (Forrest et al. 2018).
#' Set the bounds of M in the \code{M_bounds} argument, a length-2 vector where the first entry is M0, the M as B/B0 >= 1,
#' and the second entry is M1, the M as B/B0 approaches zero. Note that M0 can be greater than M1 (compensatory) 
#' or M0 can be less than M1 (depensatory).
#' 
#' @inheritParams SCA
#' @param M_bounds A numeric vector of length 2 to indicate the M as B/B0 approaches zero and one, respectively.
#' By default, set to 75\% and 125\%, respectively, of Data@@Mort[x].
#' @details
#' See \link{SCA} for more information on all arguments.
#' @references
#' Forrest, R.E., Holt, K.R., and Kronlund, A.R. 2018. Performance of alternative harvest control rules for two Pacific groundfish
#' stocks with uncertain natural mortality: Bias, robustness and trade-offs. Fisheries Research 2016: 259-286.
#'
#' @section Online Documentation:
#' Model description and equations are available on the openMSE 
#' \href{https://openmse.com/features-assessment-models/2-sca/}{website}.
#' @examples
#' res <- SCA_DDM(Data = MSEtool::SimulatedData)
#' 
#' @author Q. Huynh
#' @return An object of class \linkS4class{Assessment}.
#' @seealso \link{SCA} \link{SCA_RWM} \link{plot.Assessment} \link{summary.Assessment} \link{retrospective} \link{profile} \link{make_MP}
#' @export
SCA_DDM <- function(x = 1, Data, AddInd = "B", SR = c("BH", "Ricker", "none"), 
                    vulnerability = c("logistic", "dome"), catch_eq = c("Baranov", "Pope"),
                    CAA_dist = c("multinomial", "lognormal"),
                    CAA_multiplier = 50, rescale = "mean1", max_age = Data@MaxAge,
                    start = NULL, prior = list(), fix_h = TRUE, fix_F_equilibrium = TRUE, fix_omega = TRUE, fix_tau = TRUE,
                    LWT = list(), early_dev = c("comp_onegen", "comp", "all"), late_dev = "comp50", M_bounds = NULL,
                    integrate = FALSE, silent = TRUE, opt_hess = FALSE, n_restart = ifelse(opt_hess, 0, 1),
                    control = list(iter.max = 2e5, eval.max = 4e5), inner.control = list(), ...) {
  
  out <- SCA_(x = x, Data = Data, AddInd = AddInd, SR = SR, vulnerability = vulnerability, catch_eq = catch_eq, comp = "age",
              comp_dist = CAA_dist, comp_multiplier = CAA_multiplier, rescale = rescale, max_age = max_age,
              start = start, prior = prior, fix_h = fix_h, fix_F_equilibrium = fix_F_equilibrium, fix_omega = fix_omega, fix_tau = fix_tau, 
              tv_M = "DD", M_bounds = M_bounds,
              LWT = LWT, early_dev = early_dev, late_dev = late_dev, integrate = integrate,
              silent = silent, opt_hess = opt_hess, n_restart = n_restart,
              control = control, inner.control = inner.control, ...) 
  return(out)
}
class(SCA_DDM) <- "Assess"

