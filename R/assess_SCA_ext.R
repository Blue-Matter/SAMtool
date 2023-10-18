
#' Age-structured model using fishery length composition
#' 
#' A single-fleet assessment that fits to catch, indices of abundance, and fishery length compositions. See [SCA] for all details.
#' 
#' @inheritParams SCA
#' @param CAL_dist Character, the statistical distribution for the likelihood of the catch-at-length.
#' @param CAL_multiplier Numeric for data weighting of catch-at-length matrix if `CAL_hist = "multinomial"`. A value smaller than one
#' rescales annual sample sizes to this fraction of the original sample size. Values greater than one generates a cap of the annual sample
#' size to this value.
#' 
#' @section Online Documentation:
#' Model description and equations are available on the openMSE 
#' [website](https://openmse.com/features-assessment-models/2-sca/).
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


#' SCA models with time-varying natural mortality
#' 
#' A modification of `SCA` that incorporates density-dependent effects on M based on biomass depletion (Forrest et al. 2018).
#' Set the bounds of M in the `M_bounds` argument, a length-2 vector where the first entry is M0, the M as B/B0 >= 1,
#' and the second entry is M1, the M as B/B0 approaches zero. Note that M0 can be greater than M1 (compensatory) 
#' or M0 can be less than M1 (depensatory).
#' 
#' @inheritParams SCA
#' @param M_bounds A numeric vector of length 2 to indicate the M as B/B0 approaches zero and one, respectively.
#' By default, set to 75% and 125%, respectively, of `Data@@Mort[x]`.
#' @details
#' See [SCA] for more information on all arguments.
#' @references
#' Forrest, R.E., Holt, K.R., and Kronlund, A.R. 2018. Performance of alternative harvest control rules for two Pacific groundfish
#' stocks with uncertain natural mortality: Bias, robustness and trade-offs. Fisheries Research 2016: 259-286.
#'
#' @section Online Documentation:
#' Model description and equations are available on the openMSE 
#' [website](https://openmse.com/features-assessment-models/2-sca/).
#' @examples
#' res <- SCA_DDM(Data = MSEtool::SimulatedData)
#' 
#' @author Q. Huynh
#' @return An object of class [Assessment-class].
#' @seealso [SCA] [SCA_RWM] [plot.Assessment] [summary.Assessment] [retrospective] [profile] [make_MP]
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


#' SCA with random walk in M
#'
#' `SCA_RWM` is a modification of [SCA] that incorporates a random walk in M in logit space (constant with age). 
#' Set the variance (`start$tau_M`) to a small value (0.001) in order to fix M for all years, which is functionally equivalent to [SCA].
#' 
#' @inheritParams SCA 
#' @param refyear An expression for the year for which M is used to report MSY and unfished reference points. By default, terminal year. If multiple
#' years are provided, then the mean M over the specified time period is used.
#' @param M_bounds A numeric vector of length 2 to indicate the minimum and maximum M in the random walk as a proportion of the starting M
#' (`start$M`). The default min and max are 75% and 125%, respectively.
#' @details
#' The model estimates year-specific M (constant with age) as a random walk in logit space, bounded by
#' a proportion of `start$M` (specified in `M_bounds`).
#' 
#' The starting value for the first year M (start$M) is `Data@@Mort[x]` and is fixed, unless a prior is provided (`prior$M`). 
#' The fixed SD of the random walk (`tau_M`) is 0.05, by default. 
#' 
#' Steepness and unfished recruitment in the estimation model, along with unfished reference points, correspond to spawners per recruit using the first year M. 
#' With argument `refyear`, new unfished reference points and steepness values are calculated. See examples. 
#' 
#' Alternative values can be provided in the start list (see examples):
#' \itemize{
#' \item `R0` Unfished recruitment, except when `SR = "none"` where it is mean recruitment. 
#' By default, 150% `Data@@OM$R0[x]` is used as the start value in closed-loop simulation, and 400\% of mean catch otherwise.
#' \item `h` Steepness. Otherwise, `Data@@steep[x]` is used, or 0.9 if empty.
#' \item `M` Natural mortality in the first year. Otherwise, `Data@@Mort[x]` is used.
#' \item `vul_par` Vulnerability parameters, see next paragraph.
#' \item `F` A vector of length nyears for year-specific fishing mortality.
#' \item `F_equilibrium` Equilibrium fishing mortality leading into first year of the model (to determine initial depletion). By default, 0.
#' \item `omega` Lognormal SD of the catch (observation error) when `catch_eq = "Baranov"`. By default, `Data@@CV_Cat[x]`.
#' \item `tau` Lognormal SD of the recruitment deviations (process error). By default, `Data@@sigmaR[x]`.
#' \item `tau_M` The fixed SD of the random walk in M. By default, 0.05. 
#' }
#' 
#' See [SCA] for all other information about the structure and setup of the model.
#' 
#' The SCA builds in a stock-recruit relationship into the model. Annual unfished and MSY reference points are 
#' calculated and reported in TMB_report of the [Assessment-class] object.
#'
#' @section Online Documentation:
#' Model description and equations are available on the openMSE 
#' [website](https://openmse.com/features-assessment-models/2-sca/).
#' @examples
#' res <- SCA_RWM(Data = MSEtool::SimulatedData, start = list(M_start = 0.4, tau_M = 0.05))
#' res2 <- SCA(Data = MSEtool::SimulatedData)
#' res3 <- SCA_RWM(Data = MSEtool::SimulatedData, start = list(M_start = 0.4, tau_M = 0.001))
#' 
#' # Use mean M in most recent 5 years for reporting reference points 
#' res_5r <- SCA_RWM(Data = MSEtool::SimulatedData, 
#'                   refyear = expression(seq(length(Data@Year) - 4, length(Data@Year))),
#'                   start = list(M_start = 0.4, tau_M = 0.001))
#' res_5r@SSB0 # SSB0 reported (see also res_5r@TMB_report$new_E0)
#' res_5r@TMB_report$E0 # SSB0 of Year 1 M
#' 
#' \donttest{
#' compare_models(res, res2, res3)
#' }
#' @author Q. Huynh
#' @return An object of class [Assessment-class].
#' @seealso [SCA] [SCA_DDM]
#' @export
SCA_RWM <- function(x = 1, Data, AddInd = "B", SR = c("BH", "Ricker", "none"), 
                    vulnerability = c("logistic", "dome"), catch_eq = c("Baranov", "Pope"),
                    CAA_dist = c("multinomial", "lognormal"),
                    CAA_multiplier = 50, rescale = "mean1", max_age = Data@MaxAge,
                    start = NULL, prior = list(), fix_h = TRUE, fix_F_equilibrium = TRUE, fix_omega = TRUE, fix_tau = TRUE,
                    LWT = list(), early_dev = c("comp_onegen", "comp", "all"), late_dev = "comp50", 
                    refyear = expression(length(Data@Year)), M_bounds = NULL,
                    integrate = FALSE, silent = TRUE, opt_hess = FALSE, n_restart = ifelse(opt_hess, 0, 1),
                    control = list(iter.max = 2e5, eval.max = 4e5), inner.control = list(), ...) {
  
  out <- SCA_(x = x, Data = Data, AddInd = AddInd, SR = SR, vulnerability = vulnerability, catch_eq = catch_eq, comp = "age",
              comp_dist = CAA_dist, comp_multiplier = CAA_multiplier, rescale = rescale, max_age = max_age,
              start = start, prior = prior, fix_h = fix_h, fix_F_equilibrium = fix_F_equilibrium, fix_omega = fix_omega, fix_tau = fix_tau, 
              tv_M = "walk", M_bounds = M_bounds, refyear = refyear,
              LWT = LWT, early_dev = early_dev, late_dev = late_dev, integrate = integrate,
              silent = silent, opt_hess = opt_hess, n_restart = n_restart,
              control = control, inner.control = inner.control, ...) 
  return(out)
}
class(SCA_RWM) <- "Assess"

