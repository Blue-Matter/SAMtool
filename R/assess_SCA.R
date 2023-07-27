
#' Statistical catch-at-age (SCA) model
#'
#' A generic statistical catch-at-age model (single fleet, single season) that uses catch, index, and catch-at-age composition
#' data. \code{SCA} parameterizes R0 and steepness as leading productivity parameters in the assessment model. Recruitment is estimated
#' as deviations from the resulting stock-recruit relationship. In \code{SCA2}, the mean recruitment in the time series is estimated and
#' recruitment deviations around this mean are estimated as penalized parameters (\code{SR = "none"}, similar to Cadigan 2016). The standard deviation is set high
#' so that the recruitment is almost like free parameters. Unfished and MSY reference points are not estimated, it is recommended to use yield per recruit
#' or spawning potential ratio in harvest control rules. \code{SCA_Pope} is a variant of \code{SCA} that fixes the expected catch to the observed
#' catch, and Pope's approximation is used to calculate the annual exploitation rate (U; i.e., \code{catch_eq = "Pope"}).
#'
#' @aliases SCA2 SCA_Pope
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
#' @param prior A named list for the parameters of any priors to be added to the model. See below.
#' @param fix_h Logical, whether to fix steepness to value in \code{Data@@steep} in the model for \code{SCA}. This only affects
#' calculation of reference points for \code{SCA2}.
#' @param fix_F_equilibrium Logical, whether the equilibrium fishing mortality prior to the first year of the model
#' is estimated. If \code{TRUE}, \code{F_equilibrium} is fixed to value provided in \code{start} (if provided),
#' otherwise, equal to zero (assumes unfished conditions).
#' @param fix_omega Logical, whether the standard deviation of the catch is fixed. If \code{TRUE},
#' omega is fixed to value provided in \code{start} (if provided), otherwise, value based on \code{Data@@CV_Cat}.
#' @param fix_tau Logical, the standard deviation of the recruitment deviations is fixed. If \code{TRUE},
#' tau is fixed to value provided in \code{start} (if provided), otherwise, value based on \code{Data@@sigmaR}.
#' @param LWT A named list (Index, CAA, Catch) of likelihood weights for the data components. For the index, a vector of length survey. For
#' CAL and Catch, a single value.
#' @param early_dev Numeric or character string describing the years for which recruitment deviations are estimated in \code{SCA}. By default,
#' equal to \code{"comp_onegen"}, where rec devs are estimated one full generation prior to the first year when catch-at-age (CAA) data are available.
#' With \code{"comp"}, rec devs are estimated starting in the first year with CAA. With \code{"all"}, rec devs start at the beginning of the model.
#' If numeric, the number of years after the first year of the model for which to start estimating rec devs. Use negative numbers for years prior to the first year.
#' @param late_dev Typically, a numeric for the number of most recent years in which recruitment deviations will
#' not be estimated in \code{SCA} (recruitment in these years will be based on the mean predicted by stock-recruit relationship).
#' By default, \code{"comp50"} uses the number of ages (smaller than the mode)
#' for which the catch-at-age matrix has less than half the abundance than that at the mode.
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
#' The basic data inputs are catch (by weight), index (by weight/biomass), and catch-at-age matrix (by numbers).
#' 
#' With \code{catch_eq = "Baranov"} (default in SCA and SCA2), annual F's are estimated parameters assuming continuous fishing over the year, while
#' an annual exploitation rate from pulse fishing in the middle of the year is estimated in \code{SCA_Pope} or \code{SCA(catch_eq = "Pope")}. 
#'
#' The annual sample sizes of the catch-at-age matrix is provided to the model (used in the likelihood for catch-at-age assuming
#' a multinomial distribution) and is manipulated via argument \code{CAA_multiplier}. This argument is
#' interpreted in two different ways depending on the value provided. If \code{CAA_multiplier > 1}, then this value will cap the annual sample sizes
#' to that number. If \code{CAA_multiplier <= 1}, then all the annual samples sizes will be re-scaled by that number, e.g. \code{CAA_multiplier = 0.1}
#' multiplies the sample size to 10\% of the original number. By default, sample sizes are capped at 50.
#'
#' Alternatively, a lognormal distribution with inverse proportion variance can be used for the catch at age (Punt and Kennedy, 1994, as
#' cited by Maunder 2011).
#'
#' For \code{start} (optional), a named list of starting values of estimates can be provided for:
#' \itemize{
#' \item \code{R0} Unfished recruitment, except when \code{SR = "none"} where it is mean recruitment. 
#' By default, 150\% Data@@OM$R0[x] is used as the start value in closed-loop simulation, and 400\% of mean catch otherwise.
#' \item \code{h} Steepness. Otherwise, Data@@steep[x] is used, or 0.9 if empty.
#' \item \code{M} Natural mortality. Otherwise, Data@@Mort[x] is used.
#' \item \code{vul_par} Vulnerability parameters, see next paragraph.
#' \item \code{F} A vector of length nyears for year-specific fishing mortality.
#' \item \code{F_equilibrium} Equilibrium fishing mortality leading into first year of the model (to determine initial depletion). By default, 0.
#' \item \code{U_equilibrium} Same as F_equilibrium when \code{catch_eq = "Pope"}. By default, 0.
#' \item \code{omega} Lognormal SD of the catch (observation error) when \code{catch_eq = "Baranov"}. By default, Data@@CV_Cat[x].
#' \item \code{tau} Lognormal SD of the recruitment deviations (process error). By default, Data@@sigmaR[x].
#' }
#' 
#' Vulnerability can be specified to be either logistic or dome. If logistic, then the parameter
#' vector \code{vul_par} is of length 2:
#' \itemize{
#' \item \code{vul_par[1]} corresponds to \code{a_95}, the age of 95\% vulnerability. \code{a_95} is a transformed parameter via logit transformation to constrain \code{a_95} to less than 75\%
#' of the maximum age: \code{a_95 = 0.75 * max_age * plogis(x[1])}, where \code{x} is the estimated vector.
#' \item \code{vul_par[2]} corresponds to \code{a_50}, the age of 50\% vulnerability. Estimated as an offset, i.e., \code{a_50 = a_95 - exp(x[2])}.
#' }
#' 
#' With dome vulnerability, a double Gaussian parameterization is used, where \code{vul_par}
#' is an estimated vector of length 4:
#' \itemize{
#' \item \code{vul_par[1]} corresponds to  \code{a_asc}, the first age of full vulnerability for the ascending limb. In the model, \code{a_asc} is estimated via logit transformation
#' to constrain \code{a_95} to less than 75\% of the maximum age: \code{a_asc = 0.75 * maxage * plogis(x[1])}, where \code{x} is the estimated vector.
#' \item \code{vul_par[2]} corresponds to \code{a_50}, the age of 50\% vulnerability for the ascending limb. Estimated as an offset, i.e.,
#' \code{a_50 = a_asc - exp(x[2])}.
#' \item \code{vul_par[3]} corresponds to \code{a_des}, the last age of full vulnerability (where the descending limb starts). Generated via logit transformation
#' to constrain between \code{a_asc} and \code{max_age}, i.e., \code{a_des = (max_age - a_asc) * plogis(x[3]) + a_asc}. By default, fixed to a small value so that the dome is effectively
#' a three-parameter function.
#' \item \code{vul_par[4]} corresponds to \code{vul_max}, the vulnerability at the maximum age. Estimated in logit space: \code{vul_max = plogis(x[4])}.
#' }
#' Vague priors of \code{vul_par[1] ~ N(0, sd = 3)}, \code{vul_par[2] ~ N(0, 3)}, \code{vul_par[3] ~ Beta(1.01, 1.01)} are used to aid convergence when parameters may not be well estimated,
#' for example, when vulnerability >> 0.5 for the youngest age class.
#' 
#' @section Online Documentation:
#' Model description and equations are available on the openMSE 
#' \href{https://openmse.com/features-assessment-models/2-sca/}{website}.
#' 
#' @references
#' Cadigan, N.G. 2016. A state-space stock assessment model for northern cod, including under-reported catches and
#' variable natural mortality rates. Canadian Journal of Fisheries and Aquatic Science 72:296-308.
#'
#' Maunder, M.N. 2011. Review and evaluation of likelihood functions for composition data in
#' stock-assessment models: Estimating the effective sample size. Fisheries Research 209:311-319.
#'
#' Punt, A.E. and Kennedy, R.B. 1997. Population modelling of Tasmanian rock lobster, Jasus edwardsii, resources. Marine and Freshwater
#' Research 48:967-980.
#'
#' @examples
#' res <- SCA(Data = MSEtool::SimulatedData)
#' res2 <- SCA2(Data = MSEtool::SimulatedData)
#' 
#' # Downweight the index
#' res3 <- SCA(Data = MSEtool::SimulatedData, LWT = list(Index = 0.1, CAA = 1))
#'
#' compare_models(res, res2)
#' @section Required Data:
#' \itemize{
#' \item \code{SCA}, \code{SCA_Pope}, and \code{SCA_Pope}: Cat, Ind, Mort, L50, L95, CAA, vbK, vbLinf, vbt0, wla, wlb, MaxAge
#' }
#' @section Optional Data:
#' \itemize{
#' \item \code{SCA}: Rec, steep, sigmaR, CV_Ind, CV_Cat
#' \item \code{SC2}: Rec, steep, CV_Ind, CV_Cat
#' \item \code{SCA_Pope}: Rec, steep, sigmaR, CV_Ind
#' }
#' @author Q. Huynh
#' @return An object of class \linkS4class{Assessment}.
#' @seealso \link{plot.Assessment} \link{summary.Assessment} \link{retrospective} \link{profile} \link{make_MP}
#' @export
SCA <- function(x = 1, Data, AddInd = "B", SR = c("BH", "Ricker", "none"), 
                vulnerability = c("logistic", "dome"), catch_eq = c("Baranov", "Pope"),
                CAA_dist = c("multinomial", "lognormal"),
                CAA_multiplier = 50, rescale = "mean1", max_age = Data@MaxAge,
                start = NULL, prior = list(), fix_h = TRUE, fix_F_equilibrium = TRUE, fix_omega = TRUE, fix_tau = TRUE,
                LWT = list(), early_dev = c("comp_onegen", "comp", "all"), late_dev = "comp50", integrate = FALSE,
                silent = TRUE, opt_hess = FALSE, n_restart = ifelse(opt_hess, 0, 1),
                control = list(iter.max = 2e5, eval.max = 4e5), inner.control = list(), ...) {
  
  out <- SCA_(x = x, Data = Data, AddInd = AddInd, SR = SR, vulnerability = vulnerability, catch_eq = catch_eq, comp = "age",
              comp_dist = CAA_dist, comp_multiplier = CAA_multiplier, 
              rescale = rescale, max_age = max_age, start = start, prior = prior, fix_h = fix_h, 
              fix_F_equilibrium = fix_F_equilibrium, fix_omega = fix_omega, fix_tau = fix_tau, 
              LWT = LWT, early_dev = early_dev, late_dev = late_dev, integrate = integrate,
              silent = silent, opt_hess = opt_hess, n_restart = n_restart, control = control, 
              inner.control = inner.control, ...)
  return(out)
}
class(SCA) <- "Assess"

#' @rdname SCA
#' @param common_dev Typically, a numeric for the number of most recent years in which a common recruitment deviation will
#' be estimated (in \code{SCA2}, uninformative years will have a recruitment closer to the mean, which can be very misleading,
#' especially near the end of the time series). By default, \code{"comp50"} uses the number of ages (smaller than the mode)
#' for which the catch-at-age matrix has less than half the abundance than that at the mode.
#' @export
SCA2 <- function(x = 1, Data, AddInd = "B", vulnerability = c("logistic", "dome"), CAA_dist = c("multinomial", "lognormal"),
                 CAA_multiplier = 50, rescale = "mean1", max_age = Data@MaxAge, start = NULL, prior = list(),
                 fix_h = TRUE, fix_F_equilibrium = TRUE, fix_omega = TRUE, fix_tau = TRUE, LWT = list(),
                 common_dev = "comp50", integrate = FALSE, silent = TRUE, opt_hess = FALSE, n_restart = ifelse(opt_hess, 0, 1),
                 control = list(iter.max = 2e5, eval.max = 4e5), inner.control = list(), ...) {
  
  out <- SCA_(x = x, Data = Data, AddInd = AddInd, SR = "none", vulnerability = vulnerability, catch_eq = "Baranov", comp = "age",
              comp_dist = CAA_dist, comp_multiplier = CAA_multiplier, rescale = rescale, max_age = max_age,
              start = start, prior = prior, fix_h = fix_h, fix_F_equilibrium = fix_F_equilibrium, fix_omega = fix_omega, fix_tau = fix_tau,
              LWT = LWT, early_dev = "all", late_dev = common_dev, integrate = integrate,
              silent = silent, opt_hess = opt_hess, n_restart = n_restart, control = control, 
              inner.control = inner.control, ...)
  return(out)
}
class(SCA2) <- "Assess"


#' @rdname SCA
#' @param fix_U_equilibrium Logical, same as `fix_F_equilibrium` for `SCA_Pope`.
#' @export
SCA_Pope <- function(x = 1, Data, AddInd = "B", SR = c("BH", "Ricker", "none"), vulnerability = c("logistic", "dome"), CAA_dist = c("multinomial", "lognormal"),
                     CAA_multiplier = 50, rescale = "mean1", max_age = Data@MaxAge, start = NULL, prior = list(),
                     fix_h = TRUE, fix_U_equilibrium = TRUE, fix_tau = TRUE, LWT = list(),
                     early_dev = c("comp_onegen", "comp", "all"), late_dev = "comp50", integrate = FALSE,
                     silent = TRUE, opt_hess = FALSE, n_restart = ifelse(opt_hess, 0, 1),
                     control = list(iter.max = 2e5, eval.max = 4e5), inner.control = list(), ...) {
  
  out <- SCA_(x = x, Data = Data, AddInd = AddInd, SR = SR, vulnerability = vulnerability, catch_eq = "Pope", comp = "age",
              comp_dist = CAA_dist, comp_multiplier = CAA_multiplier, rescale = rescale, max_age = max_age,
              start = start, prior = prior, fix_h = fix_h, fix_F_equilibrium = fix_U_equilibrium, fix_omega = TRUE, fix_tau, 
              LWT = LWT, early_dev = early_dev, late_dev = late_dev, integrate = integrate,
              silent = silent, opt_hess = opt_hess, n_restart = n_restart, control = control, 
              inner.control = inner.control, ...)
  return(out)
}
class(SCA_Pope) <- "Assess"


#' @useDynLib SAMtool
SCA_ <- function(x = 1, Data, AddInd = "B", SR = c("BH", "Ricker", "none"), 
                 vulnerability = c("logistic", "dome"), catch_eq = c("Baranov", "Pope"), comp = c("age", "length"),
                 comp_dist = c("multinomial", "lognormal"), comp_multiplier = c(50, 50), rescale = "mean1", max_age = Data@MaxAge,
                 start = NULL, prior = list(), fix_h = TRUE, fix_F_equilibrium = TRUE, fix_omega = TRUE, fix_tau = TRUE,
                 tv_M = c("none", "walk", "DD"), M_bounds = NULL, refyear = 1,
                 LWT = list(), early_dev = c("comp_onegen", "comp", "all"), late_dev = "comp50", integrate = FALSE,
                 silent = TRUE, opt_hess = FALSE, n_restart = ifelse(opt_hess, 0, 1),
                 control = list(iter.max = 2e5, eval.max = 4e5), inner.control = list(), ...) {
  
  dependencies <- "Data@Cat, Data@Ind, Data@Mort, Data@L50, Data@L95, Data@CAA, Data@vbK, Data@vbLinf, Data@vbt0, Data@wla, Data@wlb, Data@MaxAge"

  dots <- list(...)
  start <- lapply(start, eval, envir = environment())
  
  max_age <- as.integer(min(max_age, Data@MaxAge))
  n_age <- max_age + 1
  vulnerability <- match.arg(vulnerability)
  catch_eq <- match.arg(catch_eq)
  comp <- match.arg(comp, several.ok = TRUE)
  comp_dist <- match.arg(comp_dist)
  if (length(comp_multiplier) == 1) comp_multiplier <- rep(comp_multiplier, 2)
  SR <- match.arg(SR)  
  tv_M <- match.arg(tv_M)
  
  if (is.character(early_dev)) early_dev <- match.arg(early_dev)
  if (is.numeric(early_dev)) stopifnot(early_dev < length(Data@Year))
  
  if (any(names(dots) == "yind")) {
    yind <- eval(dots$yind)
  } else {
    yind <- which(!is.na(Data@Cat[x, ]))[1]
    yind <- yind:length(Data@Cat[x, ])
  }
  Year <- Data@Year[yind]
  C_hist <- Data@Cat[x, yind]
  if (any(is.na(C_hist) | C_hist < 0)) warning("Error. Catch time series is not complete.")
  
  n_y <- length(C_hist)
  if (any(names(dots) == "M_at_age") && dots$M_at_age) {
    M <- Data@Misc$StockPars$M_ageArray[x, , n_y] * Data@Obs$Mbias[x]
    prior$M <- NULL
  } else {
    M <- rep(Data@Mort[x], n_age)
  }
  a <- Data@wla[x]
  b <- Data@wlb[x]
  Linf <- Data@vbLinf[x]
  K <- Data@vbK[x]
  t0 <- Data@vbt0[x]
  La <- Linf * (1 - exp(-K * (c(0:max_age) - t0)))
  SD_La <- La * Data@LenCV[x]
  Wa <- a * La ^ b
  A50 <- min(0.5 * max_age, iVB(t0, K, Linf, Data@L50[x]))
  A95 <- max(A50+0.5, iVB(t0, K, Linf, Data@L95[x]))
  mat_age <- c(0, 1/(1 + exp(-log(19) * (c(1:max_age) - A50)/(A95 - A50)))) # Age-0 is immature
  mat_age <- mat_age/max(mat_age)
  
  if (any(comp == "age")) {
    Data <- expand_comp_matrix(Data, "CAA") # Make sure dimensions of CAA match that in catch (nyears).
    CAA_hist <- Data@CAA[x, yind, 1:n_age]
    if (max_age < Data@MaxAge) CAA_hist[, n_age] <- rowSums(Data@CAA[x, yind, n_age:(Data@MaxAge+1)], na.rm = TRUE)
    
    if (all(is.na(CAA_hist))) warning("No age composition found in Data object", call. = FALSE)
  } else {
    CAA_hist <- matrix(0, n_y, n_age)
  }
  CAA_n_nominal <- rowSums(CAA_hist, na.rm = TRUE)
  if (comp_multiplier[1] <= 1) {
    CAA_n_rescale <- comp_multiplier[1] * CAA_n_nominal
  } else {
    CAA_n_rescale <- pmin(CAA_n_nominal, comp_multiplier[1])
  }
  
  if (any(comp == "length")) {
    CAL_hist <- Data@CAL[x, yind, ]
    CAL_sum <- colSums(CAL_hist, na.rm = TRUE) # Remove length bins with zeros for all years
    CAL_cdf <- cumsum(CAL_sum)
    CAL_ind <- which(CAL_sum > 0)[1]:which.max(CAL_cdf)[1]
    CAL_hist <- CAL_hist[, CAL_ind]
    if (all(is.na(CAL_hist))) warning("No length composition found in Data object", call. = FALSE)
    
    CAL_mids <- Data@CAL_mids[CAL_ind]
    CAL_bins <- Data@CAL_bins[CAL_ind]
    PLA <- generate_PLA(La, SD_La, CAL_bins, CAL_mids)
  } else {
    CAL_hist <- matrix(0, n_y, 1)
    PLA <- matrix(1, n_age, 1)
    CAL_bins <- CAL_mids <- 1
  }
  CAL_n_nominal <- rowSums(CAL_hist, na.rm = TRUE)
  if (comp_multiplier[2] <= 1) {
    CAL_n_rescale <- comp_multiplier[2] * CAL_n_nominal
  } else {
    CAL_n_rescale <- pmin(CAL_n_nominal, comp_multiplier[2])
  }
  
  if (early_dev == "all") {
    est_early_rec_dev <- rep(1, n_age - 1)
    est_rec_dev <- rep(1, n_y)
  } else if (early_dev == "comp") {
    est_early_rec_dev <- rep(0, n_age-1)
    if (any(comp == "age")) {
      est_rec_dev <- ifelse(1:n_y < which(CAA_n_nominal > 0)[1], 0, 1)
    } else {
      est_rec_dev <- ifelse(1:n_y < which(CAL_n_nominal > 0)[1], 0, 1)
    }
  } else if (early_dev == "comp_onegen") {
    if (any(comp == "age")) {
      istart <- which(CAA_n_nominal > 0)[1] - n_age
    } else {
      istart <- which(CAL_n_nominal > 0)[1] - n_age
    }
    if (istart < 0) {
      early_start <- n_age + istart
      est_early_rec_dev <- ifelse(2:n_age < early_start, 0, 1) %>% rev()
      est_rec_dev <- rep(1, n_y)
    } else {
      est_early_rec_dev <- rep(0, n_age - 1)
      est_rec_dev <- ifelse(1:n_y < istart, 0, 1)
    }
  } else if (is.numeric(early_dev)) {
    if (early_dev > 1) {
      est_early_rec_dev <- rep(0, n_age-1)
      est_rec_dev <- ifelse(1:n_y < early_dev, 0, 11)
    } else {
      istart <- early_dev - 1
      est_early_rec_dev <- rep(c(1, 0), c(istart, n_age - istart - 1))
      est_rec_dev <- rep(1, n_y)
    }
  }
  if (tv_M == "DD") est_early_rec_dev <- rep(0, n_age-1) # Temporary for now
  
  if (is.character(late_dev) && late_dev == "comp50") {
    if (any(comp == "age")) {
      comp_ldev <- colSums(CAA_hist, na.rm = TRUE)/max(colSums(CAA_hist, na.rm = TRUE))
    } else {
      comp_ldev <- colSums(CAL_hist, na.rm = TRUE)/max(colSums(CAL_hist, na.rm = TRUE))
    }
    comp_mode <- which.max(comp_ldev)[1]
    comp50_ind <- which(comp_ldev[1:comp_mode] <= 0.5)
    comp50_ind <- comp50_ind[length(comp50_ind)]
    if (all(comp != "age")) comp50_ind <- ceiling(LinInterp(La, 0:n_age, Data@CAL_mids[comp50_ind]))
    late_dev <- ifelse(is.na(comp50_ind), 0, comp50_ind)
  }
  if (is.numeric(late_dev) && late_dev > 0) {
    if (late_dev > length(est_rec_dev)) late_dev <- length(est_rec_dev)
    ind_late <- (length(est_rec_dev) - late_dev + 1):length(est_rec_dev)
    est_rec_dev[ind_late] <- ifelse(SR == "none", max(est_rec_dev[-ind_late]), 0)
  }
  
  if (rescale == "mean1") rescale <- 1/mean(C_hist)
  
  Ind <- lapply(AddInd, Assess_I_hist, Data = Data, x = x, yind = yind)
  I_hist <- vapply(Ind, getElement, numeric(n_y), "I_hist")
  if (is.null(I_hist) || all(is.na(I_hist))) stop("No indices found.", call. = FALSE)
  
  I_sd <- vapply(Ind, getElement, numeric(n_y), "I_sd") %>% pmax(0.05)
  I_units <- vapply(Ind, getElement, numeric(1), "I_units")
  
  I_vul <- vapply(AddInd, function(xx) {
    if (xx == "B") {
      return(rep(1, n_age))
    } else if (xx == "SSB") {
      return(mat_age)
    } else if (xx == "VB") {
      return(rep(0, n_age))
    } else {
      return(Data@AddIndV[x, suppressWarnings(as.numeric(xx)), 1:n_age])
    }
  }, numeric(n_age))
  nsurvey <- ncol(I_hist)
  
  if (!is.list(LWT)) {
    if (!is.null(LWT) && length(LWT) != nsurvey) stop("LWT needs to be a vector of length ", nsurvey)
    LWT <- list(Index = LWT)
    LWT$CAA <- LWT$CAL <- LWT$Catch <- 1
  } else {
    if (is.null(LWT$Index)) LWT$Index <- rep(1, nsurvey)
    if (is.null(LWT$CAA)) LWT$CAA <- 1
    if (is.null(LWT$CAL)) LWT$CAL <- 1
    if (is.null(LWT$Catch)) LWT$Catch <- 1 
  }
  
  # Generate priors
  prior <- make_prior(prior, nsurvey, ifelse(SR == "BH", 1, 2), msg = FALSE)
  
  # M_bounds
  if (is.null(M_bounds)) {
    if (tv_M == "none") {
      M_bounds <- c(0, 1e4)
    } else {
      M_bounds <- c(0.75, 1.25) * range(M)
    }
  }
  
  data <- list(model = "SCA", C_hist = C_hist, rescale = rescale, I_hist = I_hist,
               I_sd = I_sd, I_units = I_units, I_vul = I_vul, abs_I = rep(0, nsurvey), nsurvey = nsurvey, 
               CAA_hist = apply(CAA_hist, 1, tiny_comp) %>% t(), CAA_n = CAA_n_rescale, 
               CAL_hist = apply(CAL_hist, 1, tiny_comp) %>% t(), CAL_n = CAL_n_rescale,
               LWT = c(LWT$Index, LWT$CAA, LWT$CAL, LWT$Catch),
               n_y = n_y, n_age = n_age, n_bin = ncol(PLA), 
               M_data = 1,
               weight = Wa, PLA = PLA, mat = mat_age, vul_type = vulnerability,
               SR_type = SR, comp_dist = comp_dist, catch_eq = catch_eq,
               est_early_rec_dev = est_early_rec_dev, est_rec_dev = est_rec_dev, yindF = as.integer(0.5 * n_y),
               tv_M = tv_M, M_bounds = M_bounds, use_prior = prior$use_prior, prior_dist = prior$pr_matrix,
               sim_process_error = 0L)
  if (any(names(dots) == "M_at_age") && dots$M_at_age) data$M_data <- M
  if (data$n_bin == 1) data$CAL_hist <- t(data$CAL_hist)
  
  # Starting values
  params <- list()
  if (!is.null(start)) {
    if (!is.null(start$R0) && is.numeric(start$R0)) params$R0x <- log(start$R0[1] * rescale)
    if (!is.null(start$h) && is.numeric(start$h)) {
      if (SR == "BH") {
        h_start <- (start$h[1] - 0.2)/0.8
        params$transformed_h <- logit(h_start)
      } else if (SR == "Ricker") {
        params$transformed_h <- log(start$h[1] - 0.2)
      }
    }
    if (!is.null(start$M) && is.numeric(start$M)) params$log_M0 <- log(start$M)
    if (catch_eq == "Baranov" && !is.null(start$F_equilibrium) && is.numeric(start$F_equilibrium)) {
      params$F_equilibrium <- start$F_equilibrium
    }
    if (catch_eq == "Pope" && !is.null(start$U_equilibrium) && is.numeric(start$U_equilibrium)) {
      params$F_equilibrium <- start$U_equilibrium
    }
    if (!is.null(start$vul_par) && is.numeric(start$vul_par)) {
      if (start$vul_par[1] > 0.75 * max_age) stop("start$vul_par[1] needs to be less than 0.75 * Data@MaxAge (see help).")
      if (vulnerability == "logistic") {
        if (length(start$vul_par) < 2) stop("Two parameters needed for start$vul_par with logistic vulnerability (see help).")
        if (start$vul_par[1] <= start$vul_par[2]) stop("start$vul_par[1] needs to be greater than start$vul_par[2] (see help).")
        
        params$vul_par <- c(logit(start$vul_par[1]/max_age/0.75), log(start$vul_par[1] - start$vul_par[2]))
      }
      if (vulnerability == "dome") {
        if (length(start$vul_par) < 4) stop("Four parameters needed for start$vul_par with dome vulnerability (see help).")
        if (start$vul_par[1] <= start$vul_par[2]) stop("start$vul_par[1] needs to be greater than start$vul_par[2] (see help).")
        if (start$vul_par[3] <= start$vul_par[1] || start$vul_par[3] >= max_age) {
          stop("start$vul_par[3] needs to be between start$vul_par[1] and Data@MaxAge (see help).")
        }
        if (start$vul_par[4] <= 0 || start$vul_par[4] >= 1) stop("start$vul_par[4] needs to be between 0-1 (see help).")
        
        params$vul_par <- c(logit(start$vul_par[1]/max_age/0.75), log(start$vul_par[1] - start$vul_par[2]),
                            logit(1/(max_age - start$vul_par[1])), logit(start$vul_par[4]))
      }
    }
    if (!is.null(start$F) && is.numeric(start$F)) {
      Fstart <- numeric(n_y)
      Fstart_ind <- data$yindF + 1
      Fstart[Fstart_ind] <- log(start$F[Fstart_ind])
      Fstart[-Fstart_ind] <- log(start$F[-Fstart_ind]/Fstart[Fstart_ind])
      params$log_F_dev <- Fstart
    }
    
    if (!is.null(start$omega) && is.numeric(start$omega)) params$log_omega <- log(start$omega)
    if (!is.null(start[["tau"]]) && is.numeric(start[["tau"]])) params$log_tau <- log(start[["tau"]])
    if (!is.null(start[["tau_M"]]) && is.numeric(start[["tau_M"]])) params$log_tau_M <- log(start[["tau_M"]])
  }
  
  if (is.null(params$R0x)) {
    params$R0x <- ifelse(is.null(Data@OM$R0[x]), log(mean(data$C_hist)) + 4, log(1.5 * rescale * Data@OM$R0[x]))
  }
  if (is.null(params$transformed_h)) {
    h_start <- ifelse(!fix_h && is.na(Data@steep[x]), 0.9, Data@steep[x])
    if (SR == "BH") {
      h_start <- (h_start - 0.2)/0.8
      params$transformed_h <- logit(h_start)
    } else if (SR == "Ricker") {
      params$transformed_h <- log(h_start - 0.2)
    } else {
      params$transformed_h <- 0
    }
  }
  if (is.null(params$log_M0)) params$log_M0 <- log(M) %>% mean()
  if (is.null(params$logit_M_walk)) params$logit_M_walk <- rep(0, n_y)
  if (is.null(params$F_equilibrium)) params$F_equilibrium <- 0
  if (is.null(params$vul_par)) {
    if (any(comp == "age")) {
      comp_ldev <- colSums(CAA_hist, na.rm = TRUE)/max(colSums(CAA_hist, na.rm = TRUE))
    } else {
      comp_ldev <- colSums(CAL_hist, na.rm = TRUE)/max(colSums(CAL_hist, na.rm = TRUE))
    }
    comp_mode <- which.max(comp_ldev)[1]
    
    if (is.na(Data@LFC[x]) && is.na(Data@LFS[x]) || Data@LFC[x] > Linf || Data@LFS[x] > Linf) {
      
      if (all(comp == "length")) {
        comp_mode <- ceiling(LinInterp(La, 0:max_age, Data@CAL_mids[comp_mode]))
        if (!length(comp_mode)) comp_mode <- 1
      }
      if (vulnerability == "logistic") params$vul_par <- c(logit(comp_mode/max_age/0.75), log(1))
      if (vulnerability == "dome") {
        params$vul_par <- c(logit(comp_mode/max_age/0.75), log(1), logit(1/(max_age - comp_mode)), logit(0.5))
      }
    } else {
      A5 <- min(iVB(t0, K, Linf, Data@LFC[x]), comp_mode - 1)
      Afull <- min(iVB(t0, K, Linf, Data@LFS[x]), 0.5 * max_age)
      A5 <- min(A5, Afull - 0.5)
      A50_vul <- mean(c(A5, Afull))
      
      if (vulnerability == "logistic") params$vul_par <- c(logit(Afull/max_age/0.75), log(Afull - A50_vul))
      if (vulnerability == "dome") {
        params$vul_par <- c(logit(Afull/max_age/0.75), log(Afull - A50_vul), logit(0.1/(max_age - Afull)), logit(0.5))
      }
    }
  }
  if (is.na(params$vul_par[1])) params$vul_par[1] <- 1
  if (is.null(params$log_F_dev)) {
    Fstart <- numeric(n_y)
    Fstart[data$yindF + 1] <- log(0.75 * mean(M))
    params$log_F_dev <- Fstart
  }
  
  if (is.null(params$log_omega)) {
    sigmaC <- max(0.01, sdconv(1, Data@CV_Cat[x]), na.rm = TRUE)
    params$log_omega <- log(sigmaC)
  }
  if (is.null(params[["log_tau"]])) {
    tau_start <- ifelse(is.na(Data@sigmaR[x]), 0.6, Data@sigmaR[x])
    params$log_tau <- log(tau_start)
  }
  if (is.null(params[["log_tau_M"]])) params$log_tau_M <- log(0.05)
  
  params$log_early_rec_dev <- rep(0, n_age - 1)
  params$log_rec_dev <- rep(0, n_y)
  
  LH <- list(LAA = La, SD_LAA = SD_La, CAL_mids = CAL_mids, WAA = Wa, Linf = Linf, K = K, t0 = t0, a = a, b = b, A50 = A50, A95 = A95)
  info <- list(Year = Year, data = data, params = params, LH = LH, control = control,
               inner.control = inner.control)
  
  map <- list()
  if (catch_eq == "Baranov" && any(info$data$C_hist <= 0)) {
    ind <- info$data$C_hist <= 0
    info$params$log_F_dev[ind] <- -20
    map_logF <- length(params$log_F_dev)
    map_logF[ind] <- NA
    map_logF[!ind] <- 1:sum(!ind)
    map$log_F_dev <- factor(map_logF)
  } else if (catch_eq == "Pope") {
    map$log_F_dev <- factor(rep(NA, n_y))
  }
  if (fix_h && !prior$use_prior[2]) map$transformed_h <- factor(NA)
  if (!prior$use_prior[3]) map$log_M0 <- factor(NA)
  if (tv_M != "walk") map$logit_M_walk <- factor(rep(NA, n_y))
  if (fix_F_equilibrium) map$F_equilibrium <- factor(NA)
  if (fix_omega) map$log_omega <- factor(NA)
  if (fix_tau) map$log_tau <- factor(NA)
  map$log_tau_M <- factor(NA)
  if (any(!est_early_rec_dev)) map$log_early_rec_dev <- factor(ifelse(est_early_rec_dev, 1:sum(est_early_rec_dev), NA))
  if (any(!est_rec_dev)) map$log_rec_dev <- factor(ifelse(est_rec_dev, 1:sum(est_rec_dev), NA))
  if (vulnerability == "dome") map$vul_par <- factor(c(1, 2, NA, 3))
  
  random <- NULL
  if (integrate) random <- c("log_early_rec_dev", "log_rec_dev", "logit_M_walk")
  
  obj <- MakeADFun(data = info$data, parameters = info$params, hessian = TRUE,
                   map = map, random = random, DLL = "SAMtool", inner.control = inner.control, silent = silent)
  
  if (catch_eq == "Pope") {
    # Add starting values for rec-devs and increase R0 start value if U is too high (> 0.975)
    high_U <- try(obj$report(c(obj$par, obj$env$last.par[obj$env$random]))$penalty > 0, silent = TRUE)
    if (!is.character(high_U) && !is.na(high_U) && high_U) {
      Recruit <- try(Data@Rec[x, ], silent = TRUE)
      if (is.numeric(Recruit) && length(Recruit) == n_y && any(!is.na(Recruit))) {
        log_rec_dev <- log(Recruit/mean(Recruit, na.rm = TRUE))
        log_rec_dev[is.na(est_rec_dev) | is.na(log_rec_dev) | is.infinite(log_rec_dev)] <- 0
        info$params$log_rec_dev <- log_rec_dev
        
        obj <- MakeADFun(data = info$data, parameters = info$params, hessian = TRUE,
                         map = map, random = random, DLL = "SAMtool", inner.control = inner.control, silent = silent)
      }
      while (obj$par["R0x"] < 30 && obj$report(c(obj$par, obj$env$last.par[obj$env$random]))$penalty > 0) {
        obj$par["R0x"] <- obj$par["R0x"] + 1
      }
    }
  }
  
  mod <- optimize_TMB_model(obj, control, opt_hess, n_restart)
  opt <- mod[[1]]
  SD <- mod[[2]]
  report <- obj$report(obj$env$last.par.best)
  
  Yearplusone <- c(Year, max(Year) + 1)
  YearEarly <- (Year[1] - n_age + 1):(Year[1] - 1)
  YearDev <- c(YearEarly, Year)
  YearR <- c(YearDev, max(YearDev) + 1)
  R <- c(rev(report$R_early), report$R)
  
  Dev <- structure(c(rev(report$log_early_rec_dev), report$log_rec_dev), names = YearDev)
  report$dynamic_SSB0 <- SCA_dynamic_SSB0(obj) %>% 
    structure(names = Yearplusone)
  
  nll_report <- ifelse(is.character(opt), ifelse(integrate, NA, report$nll), opt$objective)
  Assessment <- new("Assessment", Model = "SCA", 
                    Name = Data@Name, conv = SD$pdHess,
                    B0 = report$B0, R0 = report$R0, N0 = report$N0,
                    SSB0 = report$E0, VB0 = report$VB0,
                    B = structure(report$B, names = Yearplusone),
                    B_B0 = structure(report$B/report$B0, names = Yearplusone),
                    SSB = structure(report$E, names = Yearplusone),
                    SSB_SSB0 = structure(report$E/report$E0, names = Yearplusone),
                    VB = structure(report$VB, names = Yearplusone),
                    VB_VB0 = structure(report$VB/report$VB0, names = Yearplusone),
                    R = structure(R, names = YearR),
                    N = structure(rowSums(report$N), names = Yearplusone),
                    N_at_age = report$N,
                    Selectivity = matrix(report$vul, nrow = length(Year), ncol = n_age, byrow = TRUE),
                    Obs_Catch = structure(C_hist, names = Year),
                    Obs_Index = structure(I_hist, dimnames = list(Year, paste0("Index_", 1:nsurvey))),
                    Obs_C_at_age = CAA_hist,
                    Catch = structure(report$Cpred, names = Year),
                    Index = structure(report$Ipred, dimnames = list(Year, paste0("Index_", 1:nsurvey))),
                    C_at_age = report$CAApred,
                    Dev = Dev, Dev_type = "log-Recruitment deviations",
                    NLL = structure(c(nll_report, report$nll_comp, report$prior, report$penalty),
                                    names = c("Total", paste0("Index_", 1:nsurvey), "CAA", "CAL", "Catch", "Dev", "M_dev", "Prior", "Penalty")),
                    info = info, obj = obj, opt = opt, SD = SD, TMB_report = report,
                    dependencies = dependencies)
  if (tv_M != "walk") Assessment@NLL <- Assessment@NLL[names(Assessment@NLL) != "M_dev"]
  if (all(comp != "age")) Assessment@NLL <- Assessment@NLL[names(Assessment@NLL) != "CAA"]
  if (all(comp != "length")) Assessment@NLL <- Assessment@NLL[names(Assessment@NLL) != "CAL"]
  if (catch_eq == "Pope") Assessment@NLL <- Assessment@NLL[names(Assessment@NLL) != "Catch"]
  if (SR != "none") Assessment@h <- report$h
  if (catch_eq == "Baranov") {
    Assessment@FMort <- structure(report$F, names = Year)
  } else {
    Assessment@U <- structure(report$U, names = Year)
  }
  
  if (Assessment@conv) {
    SE_Early <- as.list(SD, "Std. Error")$log_early_rec_dev %>% rev()
    SE_Main <- as.list(SD, "Std. Error")$log_rec_dev
    SE_Dev <- structure(c(SE_Early, SE_Main), names = YearDev)
    if (any(is.na(SE_Dev))) {
      Dev <- Dev[seq(which(!is.na(SE_Dev))[1], length(YearDev))]
      SE_Dev <- SE_Dev[seq(which(!is.na(SE_Dev))[1], length(YearDev))]
      SE_Dev[is.na(SE_Dev)] <- 0
    }
    
    refyear <- eval(refyear)
    
    ref_pt <- ref_pt_SCA(y = refyear, obj = obj, report = report)
    if (catch_eq == "Baranov") {
      report$FMSY <- ref_pt$FMSY
    } else {
      report$UMSY <- ref_pt$UMSY
    }
    report$MSY <- ref_pt$MSY
    report$VBMSY <- ref_pt$VBMSY
    report$RMSY <- ref_pt$RMSY
    report$BMSY <- ref_pt$BMSY
    report$EMSY <- ref_pt$EMSY
    report$refyear <- refyear
    
    if (!all(refyear == 1)) { # New reference points based on change in M
      report$new_B0 <- Assessment@B0 <- ref_pt$new_B0
      report$new_E0 <- Assessment@SSB0 <- ref_pt$new_E0
      report$new_VB0 <- Assessment@VB0 <- ref_pt$new_VB0
      report$new_R0 <- Assessment@R0 <- ref_pt$new_R0
      report$new_h <- Assessment@h <- ref_pt$new_h
      
      Assessment@B_B0 <- Assessment@B/Assessment@B0
      Assessment@SSB_SSB0 <- Assessment@SSB/Assessment@SSB0
      Assessment@VB_VB0 <- Assessment@VB/Assessment@VB0
    }
    
    if (catch_eq == "Baranov") {
      Assessment@FMSY <- report$FMSY
      Assessment@F_FMSY <- structure(report$F/Assessment@FMSY, names = Year)
    } else {
      Assessment@UMSY <- report$UMSY
      Assessment@U_UMSY <- structure(report$U/Assessment@UMSY, names = Year)
    }
    Assessment@MSY <- report$MSY
    Assessment@BMSY <- report$BMSY
    Assessment@SSBMSY <- report$EMSY
    Assessment@VBMSY <- report$VBMSY
    Assessment@B_BMSY <- structure(report$B/Assessment@BMSY, names = Yearplusone)
    Assessment@SSB_SSBMSY <- structure(report$E/Assessment@SSBMSY, names = Yearplusone)
    Assessment@VB_VBMSY <- structure(report$VB/Assessment@VBMSY, names = Yearplusone)
    Assessment@Dev <- Dev
    Assessment@SE_Dev <- SE_Dev
    Assessment@TMB_report <- report
    
    catch_eq_fn <- function(Ftarget) {
      projection_SCA(Assessment, Ftarget = Ftarget, p_years = 1, p_sim = 1, 
                     obs_error = list(array(1, c(1, 1, nsurvey)), matrix(1, 1, 1)),
                     process_error = matrix(1, 1, 1)) %>% slot("Catch") %>% as.vector()
    }
    Assessment@forecast <- list(per_recruit = ref_pt[["per_recruit"]], catch_eq = catch_eq_fn)
  }
  return(Assessment)
}

ref_pt_SCA <- function(y = 1, obj, report) {

  if (obj$env$data$SR_type == "none") {
    # Fit BH 
    R0_start <- log(mean(report$R))
    h_start <- logit((0.7 - 0.2)/0.8)
    opt <- nlminb(c(R0_start, h_start), get_SR, E = report$E, R = report$R, EPR0 = report$EPR0)
    SR_par <- get_SR(opt$par, E = report$E, R = report$R, EPR0 = report$EPR0, opt = FALSE)
    
    Arec <- SR_par$Arec
    Brec <- SR_par$Brec
    SR <- "BH"
  } else {
    SR_par <- NULL
    SR <- obj$env$data$SR_type
    Arec <- report$Arec
    Brec <- report$Brec
  }
  
  M <- apply(report$M[y, , drop = FALSE], 2, mean)
  weight <- obj$env$data$weight
  mat <- obj$env$data$mat
  vul <- report$vul
  
  catch_eq <- obj$env$data$catch_eq
  
  tv_M <- obj$env$data$tv_M
  M_bounds <- obj$env$data$M_bounds
  B0 <- report$B0
  
  max_F <- ifelse(catch_eq == "Baranov", 4, 0.99)
  opt2 <- optimize(yield_fn_SCA, interval = c(1e-4, max_F), M = M, mat = mat, weight = weight, vul = vul, 
                   SR = SR, Arec = Arec, Brec = Brec, catch_eq = catch_eq, B0 = B0, tv_M = tv_M, M_bounds = M_bounds)
  opt3 <- yield_fn_SCA(opt2$minimum, M = M, mat = mat, weight = weight, vul = vul, SR = SR, 
                       Arec = Arec, Brec = Brec, opt = FALSE, catch_eq = catch_eq, B0 = B0, tv_M = tv_M, M_bounds = M_bounds)
  
  if (catch_eq == "Baranov") {
    FMSY <- opt2$minimum
  } else {
    UMSY <- opt2$minimum
  }
  MSY <- -1 * opt2$objective
  VBMSY <- opt3["VB"]
  RMSY <- opt3["R"]
  BMSY <- opt3["B"]
  EMSY <- opt3["E"]
  
  if (catch_eq == "Baranov") {
    Fvec <- seq(0, 2.5 * FMSY, length.out = 100)
  } else {
    Fvec <- seq(0, 0.99, 0.01)
  }
  
  yield <- lapply(Fvec,
                  yield_fn_SCA, M = M, mat = mat, weight = weight, vul = vul, SR = SR, 
                  Arec = Arec, Brec = Brec, opt = FALSE, catch_eq = catch_eq, B0 = B0, tv_M = tv_M, M_bounds = M_bounds)
  EPR <- vapply(yield, getElement, numeric(1), "EPR")
  YPR <- vapply(yield, getElement, numeric(1), "YPR")
  
  new_B0 <- yield[[1]]["B"] # New due to change in M
  new_E0 <- yield[[1]]["E"]
  new_VB0 <- yield[[1]]["VB"]
  new_R0 <- yield[[1]]["R"]
  if (SR == "BH") {
    new_h <- Arec * EPR[1]/ (4 + Arec * EPR[1])
  } else {
    new_h <- 0.2 * (Arec * EPR[1])^0.8
  }
  
  if (catch_eq == "Baranov") {
    return(list(FMSY = FMSY, MSY = MSY, VBMSY = VBMSY, RMSY = RMSY, BMSY = BMSY, EMSY = EMSY,
                per_recruit = data.frame(FM = Fvec, SPR = EPR/EPR[1], YPR = YPR), SR_par = SR_par,
                new_B0 = new_B0, new_E0 = new_B0, new_VB0 = new_VB0, new_R0 = new_R0, new_h = new_h))
  } else {
    return(list(UMSY = UMSY, MSY = MSY, VBMSY = VBMSY, RMSY = RMSY, BMSY = BMSY, EMSY = EMSY,
                per_recruit = data.frame(U = Fvec, SPR = EPR/EPR[1], YPR = YPR), SR_par = SR_par,
                new_B0 = new_B0, new_E0 = new_B0, new_VB0 = new_VB0, new_R0 = new_R0, new_h = new_h))
  }
  
}

yield_fn_SCA <- function(x, M, mat, weight, fec = mat * weight, vul, SR = c("BH", "Ricker"), Arec, Brec, 
                         catch_eq = c("Baranov", "Pope"), opt = TRUE, x_transform = FALSE, B0 = 1,
                         tv_M = c("none", "walk", "DD"), M_bounds = NULL, spawn_time_frac = 0) {
  if (is.null(tv_M)) tv_M <- "none"
  tv_M <- match.arg(tv_M)
  
  if (tv_M != "DD") {
    yield_fn_SCA_int(x, M, mat, weight, fec, vul, SR, Arec, Brec, catch_eq, opt, x_transform, spawn_time_frac)
  } else {
    
    dep <- M_DD <- numeric(21)
    dep[1] <- 0.4
    for(i in 1:20) {
      M_DD[i] <- ifelse(dep[i] >= 1, M_bounds[1], 
                        ifelse(dep[i] <= 0, M_bounds[2], M_bounds[1] + (M_bounds[2] - M_bounds[1]) * (1 - dep[i])))
      out <- yield_fn_SCA_int(x, M = rep(M_DD[i], length(mat)), mat, weight, fec, vul, SR, Arec, Brec, catch_eq, 
                              opt = FALSE, x_transform = x_transform, spawn_time_frac = spawn_time_frac)
      if (abs(out["B"]/B0 - dep[i]) <= 1e-4) break
      dep[i+1] <- out["B"]/B0
    }
    
    if (opt) {
      return(-1 * out["Yield"])
    } else {
      return(out)
    }
  }
}

yield_fn_SCA_int <- function(x, M, mat, weight, fec = mat * weight, vul, SR = c("BH", "Ricker"), Arec, Brec, 
                             catch_eq = c("Baranov", "Pope"), opt = TRUE, x_transform = FALSE,
                             spawn_time_frac = 0) {
  SR <- match.arg(SR)
  catch_eq <- match.arg(catch_eq)
  if (catch_eq == "Baranov") {
    FMort <- ifelse(x_transform, exp(x), x)
    Z <- vul * FMort + M
    surv <- exp(-Z)
    spawn_surv <- exp(-spawn_time_frac * Z)
  } else {
    U <- ifelse(x_transform, ilogit(x), x)
    surv <- exp(-M) * (1 - vul * U)
    spawn_surv <- exp(-spawn_time_frac * M)
    if (spawn_time_frac > 0.5) spawn_surv <- spawn_surv * (1 - vul * U) # If tied, spawning goes first
  }
  n_age <- length(M)
  NPR <- calc_NPR(surv, n_age)
  EPR <- sum(NPR * spawn_surv * fec)
  if (SR == "BH") {
    Req <- (Arec * EPR - 1)/(Brec * EPR)
  } else if (SR == "Ricker") {
    Req <- log(Arec * EPR)/(Brec * EPR)
  }
  
  if (catch_eq == "Baranov") {
    CPR <- Baranov(vul, FMort, M, NPR)
  } else {
    CPR <- vul * U * NPR * exp(-0.5 * M)
  }
  YPR <- sum(CPR * weight)
  Yield <- YPR * Req
  if (opt) {
    return(-1 * Yield)
  } else {
    
    B <- Req * sum(NPR * weight)
    E <- Req * EPR
    VB <- Req * sum(NPR * vul * weight)
    
    return(c(EPR = EPR, Yield = Yield, YPR = YPR, B = B, E = E, VB = VB, R = Req))
  }
}

SCA_dynamic_SSB0 <- function(obj, par = obj$env$last.par.best, ...) {
  
  newdata <- obj$env$data
  if (obj$env$data$catch_eq == "Pope") {
    newdata$C_hist <- rep(1e-8, newdata$n_y)
    par[names(par) == "F_equilibrium"] <- 0
    
    obj2 <- MakeADFun(data = newdata, parameters = clean_tmb_parameters(obj), map = obj$env$map, 
                      random = obj$env$random, DLL = "SAMtool", silent = TRUE)
    out <- obj2$report(par)$E
  } else {
    par[names(par) == "log_F_dev"] <- log(1e-8)
    par[names(par) == "F_equilibrium"] <- 0
    out <- obj$report(par)$E
  }
  return(out)
}

get_SR <- function(pars, E, R, EPR0, opt = TRUE, figure = FALSE, type = c("BH", "Ricker"), fix_h = FALSE, h = NULL) {
  type <- match.arg(type)
  
  R0 <- exp(pars[1])
  E0 <- R0 * EPR0
  if (type == "BH") {
    if (!fix_h) h <- 0.2 + 0.8 * ilogit(pars[2])
    Arec <- 4*h/(1-h)/EPR0
    Brec <- (5*h-1)/(1-h)/E0
    
    Rpred <- Arec * E / (1 + Brec * E)
  } else if (type == "Ricker") {
    if (!fix_h) h <- 0.2 + exp(pars[2])
    Arec <- 1/EPR0 * (5*h)^1.25
    Brec <- 1.25 * log(5*h) / E0
    
    Rpred <- Arec * E * exp(-Brec * E)
  }
  sigmaR <- sqrt(sum((log(R/Rpred))^2)/length(R))
  
  if (opt){
    return(-sum(dnorm(log(R/Rpred), 0, sigmaR, log = TRUE)))
  } else {
    
    if (figure) {
      plot(E, R, ylim = c(0, max(R, R0)), xlim = c(0, max(E, E0)), xlab = "SSB", ylab = "Recruitment")
      
      E2 <- seq(0, E0, length.out = 500)
      if (type == "BH") Rpred2 <- Arec * E2 / (1 + Brec * E2)
      if (type == "Ricker") Rpred2 <- Arec * E2 * exp(-Brec * E2)
      
      lines(E2, Rpred2, col = "blue")
      abline(v = c(0.2 * E0, E0), h = c(h * R0, R0), lty = 2, col = "red")
      legend("topright", legend = c(paste0("h = ", round(h, 3)), paste0("log(R0) = ", round(log(R0), 3))), bty = "n")
    }
    
    return(list(Rpred = Rpred, R0 = R0, h = h, E0 = E0, sigmaR = sigmaR, Arec = Arec, Brec = Brec))
  }
}

generate_PLA <- function(LAA, SD_LAA, len_bins, len_mids) {
  n_age <- length(LAA)
  n_bin <- length(len_mids)
  
  PLA <- matrix(0, n_age, n_bin)
  PLA[, n_bin] <- 1 - pnorm(len_bins[n_bin], LAA, SD_LAA)
  PLA[, 1] <- pnorm(len_bins[2], LAA, SD_LAA)
  PLA[, 3:n_bin - 1] <- vapply(1:n_age, function(x) {
    pnorm(len_bins[3:n_bin], LAA[x], SD_LAA[x]) - pnorm(len_bins[3:n_bin - 1], LAA[x], SD_LAA[x])
  }, numeric(n_bin - 2)) %>% t()
  return(PLA)
}



