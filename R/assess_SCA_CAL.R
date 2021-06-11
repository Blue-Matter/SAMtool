
#' Age-structured model using fishery length composition
#' 
#' A single-fleet assessment that fits to catch, indices of abundance, and fishery length compositions. See \link{SCA} for all details.
#' 
#' @inheritParams SCA
#' @param CAL_dist Character, the statistical distribution for the likelihood of the catch-at-length.
#' @param CAL_multiplier Numeric for data weighting of catch-at-length matrix if \code{CAL_hist = "multinomial"}. A value smaller than one
#' rescales annual sample sizes to this fraction of the original sample size. Values greater than one generates a cap of the annual sample
#' size to this value.
#' 
#' @section Online Documentation:
#' Model description and equations are available on the openMSE 
#' \href{https://openmse.com/features-assessment-models/2-sca/}{website}.
#' @author Q. Huynh
#' @export
SCA_CAL <- function(x = 1, Data, AddInd = "B", SR = c("BH", "Ricker", "none"), 
                    vulnerability = c("logistic", "dome"), catch_eq = c("Baranov", "Pope"),
                    CAL_dist = c("multinomial", "lognormal"),
                    CAL_multiplier = 50, rescale = "mean1", max_age = Data@MaxAge,
                    start = NULL, prior = list(), fix_h = TRUE, fix_F_equilibrium = TRUE, fix_omega = TRUE, fix_tau = TRUE,
                    LWT = list(), early_dev = c("comp_onegen", "comp", "all"), late_dev = "comp50", integrate = FALSE,
                    silent = TRUE, opt_hess = FALSE, n_restart = ifelse(opt_hess, 0, 1),
                    control = list(iter.max = 2e5, eval.max = 4e5), inner.control = list(), ...) {
  
  out <- SCA_(x = x, Data = Data, AddInd = AddInd, SR = SR, vulnerability = vulnerability, catch_eq = catch_eq, comp = "length",
              comp_dist = CAL_dist, comp_multiplier = CAL_multiplier, 
              rescale = rescale, max_age = max_age, start = start, prior = prior, fix_h = fix_h, 
              fix_F_equilibrium = fix_F_equilibrium, fix_omega = fix_omega, fix_tau = fix_tau, 
              LWT = LWT, early_dev = early_dev, late_dev = late_dev, integrate = integrate,
              silent = silent, opt_hess = opt_hess, n_restart = n_restart, control = control, 
              inner.control = inner.control, ...)
  return(out)
}
class(SCA_CAL) <- "Assessment"