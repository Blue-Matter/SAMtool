
#' Statistical catch-at-age (SCA) model
#'
#' A generic statistical catch-at-age model (single fleet, single season) that uses catch, index, and catch-at-age composition
#' data. \code{SCA} parameterizes R0 and steepness as leading productivity parameters in the assessment model. Recruitment is estimated
#' as deviations from the resulting stock-recruit relationship. In \code{SCA2}, the mean recruitment in the time series is estimated and
#' recruitment deviations around this mean are estimated as penalized parameters (similar to Cadigan 2016). The standard deviation is set high
#' so that the recruitment is almost like free parameters. Unfished and MSY reference points are inferred afterwards from the assessment output
#' (SSB and recruitment estimates). \code{SCA_Pope} is a variant of \code{SCA} that fixes the expected catch to the observed
#' catch, and Pope's approximation is used to calculate the annual harvest rate (U).
#'
#' @aliases SCA2 SCA_Pope
#' @param x A position in the Data object (by default, equal to one for assessments).
#' @param Data An object of class Data
#' @param SR Stock-recruit function (either \code{"BH"} for Beverton-Holt or \code{"Ricker"}).
#' @param vulnerability Whether estimated vulnerability is \code{"logistic"} or \code{"dome"} (double-normal).
#' See details for parameterization.
#' @param CAA_dist Whether a multinomial or lognormal distribution is used for likelihood of the catch-at-age matrix. See details.
#' @param CAA_multiplier Numeric for data weighting of catch-at-age matrix if \code{CAA_hist = "multinomial"}. Otherwise ignored. See details.
#' @param I_type Whether the index surveys population biomass (B; this is the default in the DLMtool operating model),
#' vulnerable biomass (VB), or spawning stock biomass (SSB).
#' @param rescale A multiplicative factor that rescales the catch in the assessment model, which
#' can improve convergence. By default, \code{"mean1"} scales the catch so that time series mean is 1, otherwise a numeric.
#' Output is re-converted back to original units.
#' @param max_age Integer, the maximum age (plus-group) in the model.
#' @param start Optional list of starting values. Entries can be expressions that are evaluated in the function. See details.
#' @param fix_h Logical, whether to fix steepness to value in \code{Data@@steep} in the model for \code{SCA}. This only affects
#' calculation of reference points for \code{SCA2}.
#' @param fix_F_equilibrium Logical, whether the equilibrium fishing mortality prior to the first year of the model
#' is estimated. If \code{TRUE}, \code{F_equilibrium} is fixed to value provided in \code{start} (if provided),
#' otherwise, equal to zero (assumes unfished conditions).
#' @param fix_U_equilibrium Logical, same as `fix_F_equilibrium` for `SCA_Pope`.
#' @param fix_omega Logical, whether the standard deviation of the catch is fixed. If \code{TRUE},
#' sigma is fixed to value provided in \code{start} (if provided), otherwise, value based on \code{Data@@CV_Cat}.
#' @param fix_sigma Logical, whether the standard deviation of the index is fixed. If \code{TRUE},
#' sigma is fixed to value provided in \code{start} (if provided), otherwise, value based on \code{Data@@CV_Ind}.
#' @param fix_tau Logical, the standard deviation of the recruitment deviations is fixed. If \code{TRUE},
#' tau is fixed to value provided in \code{start} (if provided), otherwise, value based on \code{Data@@sigmaR}.
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
#' @param integrate Logical, whether the likelihood of the model integrates over the likelihood
#' of the recruitment deviations (thus, treating it as a random effects/state-space variable).
#' Otherwise, recruitment deviations are penalized parameters.
#' @param silent Logical, passed to \code{\link[TMB]{MakeADFun}}, whether TMB
#' will print trace information during optimization. Used for dignostics for model convergence.
#' @param opt_hess Logical, whether the hessian function will be passed to \code{\link[stats]{nlminb}} during optimization
#' (this generally reduces the number of iterations to convergence, but is memory and time intensive and does not guarantee an increase
#' in convergence rate). Ignored if \code{integrate = TRUE}.
#' @param n_restart The number of restarts (calls to \code{\link[stats]{nlminb}}) in the optimization procedure, so long as the model
#' hasn't converged. The optimization continues from the parameters from the previous (re)start.
#' @param control A named list of agruments for optimization to be passed to
#' \code{\link[stats]{nlminb}}.
#' @param inner.control A named list of arguments for optimization of the random effects, which
#' is passed on to \code{\link[TMB]{newton}}.
#' @param ... Other arguments to be passed.
#' @details
#' The basic data inputs are catch (by weight), index (by weight/biomass), and catch-at-age matrix (by numbers).
#' Annual F's are estimated parameters assuming continuous fishing over the year. Note: prior to version 1.2, catches were assumed
#' to be known perfectly with an annual harvest rate from pulse fishing in \code{SCA}. That feature has now moved to \code{SCA_Pope}.
#'
#' By default, steepness is fixed in the model to the value in \code{Data@@steep}.
#'
#' The annual sample sizes of the catch-at-age matrix is provided to the model (used in the likelihood for catch-at-age assuming
#' a multinomial distribution),
#' and is manipulated via argument \code{CAA_multiplier}. This argument is
#' interpreted in two different ways depending on the value provided.
#' If \code{CAA_multiplier > 1}, then this value will cap the annual sample sizes
#' to that number. If \code{CAA_multiplier <= 1}, then all the annual samples sizes
#' will be re-scaled by that number. By default, sample sizes are capped at 50.
#'
#' Alternatively, a lognormal distribution with inverse proportion variance can be used for the catch at age (Punt and Kennedy, 1994, as
#' cited by Maunder 2011).
#'
#' For \code{start} (optional), a named list of starting values of estimates can be provided for:
#' \itemize{
#' \item \code{R0} Virgin recruitment, only for \code{SCA}.
#' \item \code{h} Steepness, only for \code{SCA}. If not provided, the value in \code{Data@@steep} is used.
#' \item \code{meanR} Mean recruitment, only for \code{SCA2}.
#' \item \code{F_equilibrium} Fishing mortality prior to the first year of model, e.g. zero means unfished conditions. Defaults to zero.
#' \item \code{vul_par} Vulnerability parameters (length 2 vector for logistic or length 4 for dome, see below). Users should provide
#' estimates of the parameters in normal space, e.g. \code{vul_max} between 0-1, and the function will perform the appropriate transformations for the model.
#' \item \code{F} A vector of F's of length nyears, \code{length(Data@@Year)}. If not provided, defaults to 0.1.
#' \item \code{omega} Standard deviation of catch. If not provided, the value based on \code{Data@@CV_Cat} is used.
#' \item \code{sigma} Standard deviation of index. If not provided, the value based on \code{Data@@CV_Ind} is used.
#' \item \code{tau} Standard deviation of recruitment deviations. If not provided, the value in \code{Data@@sigmaR} is used.
#' }
#'
#' Vulnerability can be specified to be either logistic or dome. If logistic, then the parameter
#' vector \code{vul_par} is of length 2:
#' \itemize{
#' \item \code{vul_par[1]}: \code{a_95}, the age of 95\% vulnerability, via logit transformation to constrain \code{a_95} to less than 75\%
#' of the maximum age: \code{a_95 = 0.75 * max_age * plogis(vul_par[1])}.
#' \item \code{vul_par[2]}: \code{a_50}, the age of 50\% vulnerability as an offset, i.e., \code{a_50 = a_95 - exp(vul_par[2])}.
#' }
#' A vague prior for \code{vul_par[2] ~ N(0, sd = 3)} is used to aid convergence, for example, when vulnerability >> 0.5 for the youngest age class.
#'
#' With dome vulnerability, a double Gaussian parameterization is used, where \code{vul_par}
#' is an estimated vector of length 4:
#' \itemize{
#' \item \code{vul_par[1]}: \code{a_asc}, the first age of full vulnerability for the ascending limb, via logit transformation
#' to constrain \code{a_95} to less than 75\% of the maximum age: \code{a_asc = 0.75 * maxage * plogis(vul_par[1])}.
#' \item \code{vul_par[2]}: \code{a_50}, the age of 50\% vulnerability for the ascending limb as an offset, i.e.,
#' \code{a_50 = a_asc - exp(vul_par[2])}.
#' \item \code{vul_par[3]}: \code{a_des}, the last age of full vulnerability (where the descending limb starts) via logit transformation
#' to constrain between \code{a_asc} and \code{max_age},
#' i.e., \code{a_des = (max_age - a_asc) * plogis(vul_par[3]) + a_asc}. By default, fixed to a small value so that the dome is effectively
#' a three-parameter function.
#' \item \code{vul_par[4]}: \code{vul_max}, the vulnerability (in logit space) at the maximum age.
#' }
#' Vague priors of \code{vul_par[2] ~ N(0, sd = 3)} and \code{vul_par[3] ~ N(0, 3)} are used to aid convergence,
#' for example, when vulnerability >> 0.5 for the youngest age class.
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
#' res <- SCA(Data = DLMtool::SimulatedData)
#' res2 <- SCA2(Data = DLMtool::SimulatedData)
#'
#' compare_models(res, res2)
#'
#' SCA_assess <- SCA2(Data = DLMtool::Simulation_1)
#'
#' \dontrun{
#' plot(res)
#' }
#'
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
#' @useDynLib MSEtool
#' @export
SCA <- function(x = 1, Data, SR = c("BH", "Ricker"), vulnerability = c("logistic", "dome"), CAA_dist = c("multinomial", "lognormal"),
                CAA_multiplier = 50, I_type = c("B", "VB", "SSB"), rescale = "mean1", max_age = Data@MaxAge,
                start = NULL, fix_h = TRUE, fix_F_equilibrium = TRUE, fix_omega = TRUE, fix_sigma = FALSE, fix_tau = TRUE,
                early_dev = c("comp_onegen", "comp", "all"), late_dev = "comp50", integrate = FALSE,
                silent = TRUE, opt_hess = FALSE, n_restart = ifelse(opt_hess, 0, 1),
                control = list(iter.max = 2e5, eval.max = 4e5), inner.control = list(), ...) {
  dependencies <- "Data@Cat, Data@Ind, Data@Mort, Data@L50, Data@L95, Data@CAA, Data@vbK, Data@vbLinf, Data@vbt0, Data@wla, Data@wlb, Data@MaxAge"
  dots <- list(...)
  start <- lapply(start, eval, envir = environment())

  max_age <- as.integer(min(max_age, Data@MaxAge))
  vulnerability <- match.arg(vulnerability)
  CAA_dist <- match.arg(CAA_dist)
  SR <- match.arg(SR)
  I_type <- match.arg(I_type)
  early_dev <- match.arg(early_dev)

  if(any(names(dots) == "yind")) {
    yind <- eval(dots$yind)
  } else {
    yind <- which(!is.na(Data@Cat[x, ]))[1]
    yind <- yind:length(Data@Cat[x, ])
  }
  Year <- Data@Year[yind]
  C_hist <- Data@Cat[x, yind]
  if(any(is.na(C_hist) | C_hist < 0)) warning("Error. Catch time series is not complete.")
  I_hist <- Data@Ind[x, yind]
  Data <- expand_comp_matrix(Data, "CAA") # Make sure dimensions of CAA match that in catch (nyears).
  CAA_hist <- Data@CAA[x, yind, 1:max_age]
  if(max_age < Data@MaxAge) CAA_hist[, max_age] <- rowSums(Data@CAA[x, yind, max_age:Data@MaxAge], na.rm = TRUE)

  CAA_n_nominal <- rowSums(CAA_hist)
  if(CAA_multiplier <= 1) {
    CAA_n_rescale <- CAA_multiplier * CAA_n_nominal
  } else CAA_n_rescale <- pmin(CAA_multiplier, CAA_n_nominal)

  n_y <- length(C_hist)
  M <- rep(Data@Mort[x], max_age)
  a <- Data@wla[x]
  b <- Data@wlb[x]
  Linf <- Data@vbLinf[x]
  K <- Data@vbK[x]
  t0 <- Data@vbt0[x]
  La <- Linf * (1 - exp(-K * (c(1:max_age) - t0)))
  Wa <- a * La ^ b
  A50 <- min(0.5 * max_age, iVB(t0, K, Linf, Data@L50[x]))
  A95 <- max(A50+0.5, iVB(t0, K, Linf, Data@L95[x]))
  mat_age <- 1/(1 + exp(-log(19) * (c(1:max_age) - A50)/(A95 - A50)))
  mat_age <- mat_age/max(mat_age)
  LH <- list(LAA = La, WAA = Wa, Linf = Linf, K = K, t0 = t0, a = a, b = b, A50 = A50, A95 = A95)

  if(early_dev == "all") {
    est_early_rec_dev <- rep(1, max_age-1)
    est_rec_dev <- rep(1, n_y)
  }
  if(early_dev == "comp") {
    est_early_rec_dev <- rep(0, max_age-1)
    ind1 <- which(!is.na(CAA_n_nominal))[1]
    est_rec_dev <- ifelse(1:n_y < ind1, 0, 1)
  }
  if(early_dev == "comp_onegen") {
    ind1 <- which(!is.na(CAA_n_nominal))[1] - max_age
    if(ind1 < 0) {
      early_start <- max_age + ind1
      est_early_rec_dev <- rev(ifelse(c(1:(max_age-1)) < early_start, 0, 1))
      est_rec_dev <- rep(1, n_y)
    } else {
      est_early_rec_dev <- rep(0, max_age-1)
      est_rec_dev <- ifelse(1:n_y < ind1, 0, 1)
    }
  }
  if(is.numeric(early_dev)) {
    if(early_dev > 1) {
      est_early_rec_dev <- rep(0, max_age-1)
      est_rec_dev <- ifelse(1:n_y >= early_dev, 1, 0)
    } else {
      ind1 <- early_dev - 1
      est_early_rec_dev <- c(rep(1, ind1), rep(NA, max_age-ind1-1))
      est_rec_dev <- rep(1, n_y)
    }
  }
  if(is.character(late_dev) && late_dev == "comp50") {
    CAA_all <- colSums(CAA_hist, na.rm = TRUE)/max(colSums(CAA_hist, na.rm = TRUE))
    CAA_mode <- which.max(CAA_all)[1]
    comp50_ind <- which(CAA_all[1:CAA_mode] <= 0.5)
    comp50_ind <- comp50_ind[length(comp50_ind)]
    late_dev <- ifelse(is.na(comp50_ind), 0, comp50_ind)
  }
  if(is.numeric(late_dev) && late_dev > 0) {
    if(late_dev > length(est_rec_dev)) late_dev <- length(est_rec_dev)
    ind_late <- (length(est_rec_dev) - late_dev + 1):length(est_rec_dev)
    est_rec_dev[ind_late] <- 0
  }

  if(rescale == "mean1") rescale <- 1/mean(C_hist)
  data <- list(model = "SCA", C_hist = C_hist, rescale = rescale, I_hist = I_hist,
               CAA_hist = t(apply(CAA_hist, 1, function(x) x/sum(x))),
               CAA_n = CAA_n_rescale, n_y = n_y, max_age = max_age, M = M,
               weight = Wa, mat = mat_age, vul_type = vulnerability, I_type = I_type,
               SR_type = SR, CAA_dist = CAA_dist, est_early_rec_dev = est_early_rec_dev,
               est_rec_dev = est_rec_dev, yindF = as.integer(0.5 * n_y))
  data$CAA_hist[data$CAA_hist < 1e-8] <- 1e-8

  # Starting values
  params <- list()
  if(!is.null(start)) {
    if(!is.null(start$R0) && is.numeric(start$R0)) params$R0x <- log(start$R0[1] * rescale)
    if(!is.null(start$h) && is.numeric(start$h)) {
      if(SR == "BH") {
        h_start <- (start$h[1] - 0.2)/0.8
        params$transformed_h <- logit(h_start)
      } else {
        params$transformed_h <- log(start$h[1] - 0.2)
      }
    }
    if(!is.null(start$F_equilibrium) && is.numeric(start$F_equilibrium)) params$F_equilibrium <- start$F_equilibrium
    if(!is.null(start$vul_par) && is.numeric(start$vul_par)) {
      if(start$vul_par[1] > 0.75 * max_age) stop("start$vul_par[1] needs to be less than 0.75 * Data@MaxAge (see help).")
      if(vulnerability == "logistic") {
        if(length(start$vul_par) < 2) stop("Two parameters needed for start$vul_par with logistic vulnerability (see help).")
        if(start$vul_par[1] <= start$vul_par[2]) stop("start$vul_par[1] needs to be greater than start$vul_par[2] (see help).")

        params$vul_par <- c(logit(start$vul_par[1]/max_age/0.75), log(start$vul_par[1] - start$vul_par[2]))
      }
      if(vulnerability == "dome") {
        if(length(start$vul_par) < 4) stop("Four parameters needed for start$vul_par with dome vulnerability (see help).")
        if(start$vul_par[1] <= start$vul_par[2]) stop("start$vul_par[1] needs to be greater than start$vul_par[2] (see help).")
        if(start$vul_par[3] <= start$vul_par[1] || start$vul_par[3] >= max_age) {
          stop("start$vul_par[3] needs to be between start$vul_par[1] and Data@MaxAge (see help).")
        }
        if(start$vul_par[4] <= 0 || start$vul_par[4] >= 1) stop("start$vul_par[4] needs to be between 0-1 (see help).")

        params$vul_par <- c(logit(start$vul_par[1]/max_age/0.75), log(start$vul_par[1] - start$vul_par[2]),
                            logit(1/(max_age - start$vul_par[1])), logit(start$vul_par[4]))
      }
    }
    if(!is.null(start$F) && is.numeric(start$F)) {
      Fstart <- numeric(n_y)
      Fstart_ind <- data$yindF + 1
      Fstart[Fstart_ind] <- log(start$F[Fstart_ind])
      Fstart[-Fstart_ind] <- log(start$F[-Fstart_ind]/Fstart[Fstart_ind])
      params$logF <- Fstart
    }

    if(!is.null(start$omega) && is.numeric(start$omega)) params$log_omega <- log(start$omega)
    if(!is.null(start$sigma) && is.numeric(start$sigma)) params$log_sigma <- log(start$sigma)
    if(!is.null(start$tau) && is.numeric(start$tau)) params$log_tau <- log(start$tau)
  }

  if(is.null(params$R0x)) {
    params$R0x <- ifelse(is.null(Data@OM$R0[x]), log(mean(data$C_hist)) + 4, log(1.5 * rescale * Data@OM$R0[x]))
  }
  if(is.null(params$transformed_h)) {
    h_start <- ifelse(!fix_h && is.na(Data@steep[x]), 0.9, Data@steep[x])
    if(SR == "BH") {
      h_start <- (h_start - 0.2)/0.8
      params$transformed_h <- logit(h_start)
    } else {
      params$transformed_h <- log(h_start - 0.2)
    }
  }
  if(is.null(params$F_equilibrium)) params$F_equilibrium <- 0
  if(is.null(params$vul_par)) {
    CAA_mode <- which.max(colSums(CAA_hist, na.rm = TRUE))
    if((is.na(Data@LFC[x]) && is.na(Data@LFS[x])) || (Data@LFC[x] > Linf) || (Data@LFS[x] > Linf)) {
      if(vulnerability == "logistic") params$vul_par <- c(logit(CAA_mode/max_age/0.75), log(1))
      if(vulnerability == "dome") {
        params$vul_par <- c(logit(CAA_mode/max_age/0.75), log(1), logit(1/(max_age - CAA_mode)), logit(0.5))
      }
    } else {
      A5 <- min(iVB(t0, K, Linf, Data@LFC[x]), CAA_mode-1)
      Afull <- min(iVB(t0, K, Linf, Data@LFS[x]), 0.5 * max_age)
      A5 <- min(A5, Afull - 0.5)
      A50_vul <- mean(c(A5, Afull))

      if(vulnerability == "logistic") params$vul_par <- c(logit(Afull/max_age/0.75), log(Afull - A50_vul))
      if(vulnerability == "dome") {
        params$vul_par <- c(logit(Afull/max_age/0.75), log(Afull - A50_vul), logit(0.1/(max_age - Afull)), logit(0.5))
      }
    }
  }
  if(is.na(params$vul_par[1])) params$vul_par[1] <- 1
  if(is.null(params$logF)) {
    Fstart <- numeric(n_y)
    Fstart[data$yindF + 1] <- log(0.75 * mean(data$M))
    params$logF <- Fstart
  }

  if(is.null(params$log_omega)) {
    sigmaC <- max(0.01, sdconv(1, Data@CV_Cat[x]), na.rm = TRUE)
    params$log_omega <- log(sigmaC)
  }
  if(is.null(params$log_sigma)) {
    sigmaI <- max(0.01, sdconv(1, Data@CV_Ind[x]), na.rm = TRUE)
    params$log_sigma <- log(sigmaI)
  }
  if(is.null(params$log_tau)) {
    tau_start <- ifelse(is.na(Data@sigmaR[x]), 0.6, Data@sigmaR[x])
    params$log_tau <- log(tau_start)
  }
  params$log_early_rec_dev <- rep(0, max_age - 1)
  params$log_rec_dev <- rep(0, n_y)

  info <- list(Year = Year, data = data, params = params, LH = LH, control = control,
               inner.control = inner.control)

  map <- list()
  if(any(info$data$C_hist <= 0)) {
    ind <- info$data$C_hist <= 0
    info$params$logF[ind] <- -20
    map_logF <- length(params$logF)
    map_logF[ind] <- NA
    map_logF[!ind] <- 1:sum(!ind)
    map$logF <- factor(map_logF)
  }
  if(fix_h) map$transformed_h <- factor(NA)
  if(fix_F_equilibrium) map$F_equilibrium <- factor(NA)
  if(fix_omega) map$log_omega <- factor(NA)
  if(fix_sigma) map$log_sigma <- factor(NA)
  if(fix_tau) map$log_tau <- factor(NA)
  if(any(!est_early_rec_dev)) map$log_early_rec_dev <- factor(ifelse(est_early_rec_dev, 1:sum(est_early_rec_dev), NA))
  if(any(!est_rec_dev)) map$log_rec_dev <- factor(ifelse(est_rec_dev, 1:sum(est_rec_dev), NA))
  if(vulnerability == "dome") map$vul_par <- factor(c(1, 2, NA, 3))

  random <- NULL
  if(integrate) random <- c("log_early_rec_dev", "log_rec_dev")

  obj <- MakeADFun(data = info$data, parameters = info$params, hessian = TRUE,
                   map = map, random = random, DLL = "MSEtool", inner.control = inner.control, silent = silent)

  mod <- optimize_TMB_model(obj, control, opt_hess, n_restart)
  opt <- mod[[1]]
  SD <- mod[[2]]
  report <- obj$report(obj$env$last.par.best)

  Yearplusone <- c(Year, max(Year) + 1)
  YearEarly <- (Year[1] - max_age + 1):(Year[1] - 1)
  YearDev <- c(YearEarly, Year)
  YearR <- c(YearDev, max(YearDev) + 1)
  R <- c(rev(report$R_early), report$R)

  Dev <- structure(c(rev(report$log_early_rec_dev), report$log_rec_dev), names = YearDev)

  nll_report <- ifelse(is.character(opt), ifelse(integrate, NA, report$nll), opt$objective)
  Assessment <- new("Assessment", Model = "SCA", Name = Data@Name, conv = !is.character(SD) && SD$pdHess,
                    B0 = report$B0, R0 = report$R0, N0 = report$N0,
                    SSB0 = report$E0, VB0 = report$VB0,
                    h = report$h, FMort = structure(report$F, names = Year),
                    B = structure(report$B, names = Yearplusone),
                    B_B0 = structure(report$B/report$B0, names = Yearplusone),
                    SSB = structure(report$E, names = Yearplusone),
                    SSB_SSB0 = structure(report$E/report$E0, names = Yearplusone),
                    VB = structure(report$VB, names = Yearplusone),
                    VB_VB0 = structure(report$VB/report$VB0, names = Yearplusone),
                    R = structure(R, names = YearR),
                    N = structure(rowSums(report$N), names = Yearplusone),
                    N_at_age = report$N,
                    Selectivity = matrix(report$vul, nrow = length(Year),
                                         ncol = max_age, byrow = TRUE),
                    Obs_Catch = structure(C_hist, names = Year),
                    Obs_Index = structure(I_hist, names = Year),
                    Obs_C_at_age = CAA_hist,
                    Catch = structure(report$Cpred, names = Year),
                    Index = structure(report$Ipred, names = Year),
                    C_at_age = report$CAApred,
                    Dev = Dev, Dev_type = "log-Recruitment deviations",
                    NLL = structure(c(nll_report, report$nll_comp, report$penalty + report$prior),
                                    names = c("Total", "Index", "CAA", "Catch", "Dev", "Penalty")),
                    info = info, obj = obj, opt = opt, SD = SD, TMB_report = report,
                    dependencies = dependencies)

  if(Assessment@conv) {
    ref_pt <- SCA_MSY_calc(Arec = report$Arec, Brec = report$Brec, M = M, weight = Wa, mat = mat_age, vul = report$vul, SR = SR)
    report <- c(report, ref_pt)

    if(integrate) {
      SE_Early <- ifelse(est_early_rec_dev, sqrt(SD$diag.cov.random[names(SD$par.random) == "log_early_rec_dev"]), NA)
      SE_Main <- ifelse(est_rec_dev, sqrt(SD$diag.cov.random[names(SD$par.random) == "log_rec_dev"]), NA)
    } else {
      SE_Early <- ifelse(est_early_rec_dev, sqrt(diag(SD$cov.fixed)[names(SD$par.fixed) == "log_early_rec_dev"]), NA)
      SE_Main <- ifelse(est_rec_dev, sqrt(diag(SD$cov.fixed)[names(SD$par.fixed) == "log_rec_dev"]), NA)
    }

    SE_Dev <- structure(c(rev(SE_Early), SE_Main), names = YearDev)

    first_non_zero <- which(!is.na(SE_Dev))[1]
    if(!is.na(first_non_zero) && first_non_zero > 1) {
      Dev <- Dev[-c(1:(first_non_zero - 1))]
      SE_Dev <- SE_Dev[-c(1:(first_non_zero - 1))]
      SE_Dev[is.na(SE_Dev)] <- 0
    }

    Assessment@FMSY <- report$FMSY
    Assessment@MSY <- report$MSY
    Assessment@BMSY <- report$BMSY
    Assessment@SSBMSY <- report$EMSY
    Assessment@VBMSY <- report$VBMSY
    Assessment@F_FMSY <- structure(report$F/report$FMSY, names = Year)
    Assessment@B_BMSY <- structure(report$B/report$BMSY, names = Yearplusone)
    Assessment@SSB_SSBMSY <- structure(report$E/report$EMSY, names = Yearplusone)
    Assessment@VB_VBMSY <- structure(report$VB/report$VBMSY, names = Yearplusone)
    Assessment@Dev <- Dev
    Assessment@SE_Dev <- SE_Dev
    Assessment@TMB_report <- report
  }
  return(Assessment)
}
class(SCA) <- "Assess"



SCA_MSY_calc <- function(Arec, Brec, M, weight, mat, vul, SR = c("BH", "Ricker")) {
  SR <- match.arg(SR)
  maxage <- length(M)

  solveMSY <- function(logF) {
    Fmort <- exp(logF)
    surv <- exp(-vul * Fmort - M)
    NPR <- c(1, cumprod(surv[1:(maxage-1)]))
    NPR[maxage] <- NPR[maxage]/(1 - surv[maxage])
    EPR <- sum(NPR * mat * weight)
    if(SR == "BH") Req <- (Arec * EPR - 1)/(Brec * EPR)
    if(SR == "Ricker") Req <- log(Arec * EPR)/(Brec * EPR)
    CPR <- vul * Fmort/(vul * Fmort + M) * NPR * (1 - exp(-vul * Fmort - M))
    Yield <- Req * sum(CPR * weight)
    return(-1 * Yield)
  }

  opt2 <- optimize(solveMSY, interval = c(-6, 6))
  FMSY <- exp(opt2$minimum)
  MSY <- -1 * opt2$objective

  surv_MSY <- exp(-vul * FMSY - M)
  NPR_MSY <- c(1, cumprod(surv_MSY[1:(maxage-1)]))
  NPR_MSY[maxage] <- NPR_MSY[maxage]/(1 - surv_MSY[maxage])

  EPR_MSY <- sum(NPR_MSY * weight * mat)
  if(SR == "BH") RMSY <- (Arec * EPR_MSY - 1)/(Brec * EPR_MSY)
  if(SR == "Ricker") RMSY <- log(Arec * EPR_MSY)/(Brec * EPR_MSY)

  VBMSY <- RMSY * sum(NPR_MSY * weight * vul)
  BMSY <- RMSY * sum(NPR_MSY * weight)
  EMSY <- RMSY * EPR_MSY

  return(list(FMSY = FMSY, MSY = MSY, VBMSY = VBMSY, RMSY = RMSY, BMSY = BMSY, EMSY = EMSY))
}


SCA_refpt_calc <- function(E, R, weight, mat, M, vul, SR, fix_h, h) {
  # Per-recruit quantities
  maxage <- length(M)
  surv0 <- exp(-M)
  NPR0 <- c(1, cumprod(surv0[1:(maxage-1)]))
  NPR0[maxage] <- NPR0[maxage]/(1 - exp(-M[maxage]))
  EPR0 <- sum(NPR0 * weight * mat)

  # Fit stock-recruit curve
  Rpred <- sigmaR <- NULL
  solve_SR_par <- function(x, h = NULL) {
    R0 <- exp(x[1])
    E0 <- R0 * EPR0
    if(!fix_h) {
      if(SR == "BH") h <- 0.2 + 0.8 * ilogit(x[2])
      if(SR == "Ricker") h <- 0.2 + exp(x[2])
    }
    if(SR == "BH") Rpred <<- (0.8 * R0 * h * E)/(0.2 * EPR0 * R0 *(1-h)+(h-0.2)*E)
    if(SR == "Ricker") Rpred <<- E/EPR0 * (5*h)^(1.25 * (1 - E/E0))
    sigmaR <<- sqrt(sum((log(R/Rpred))^2)/length(Rpred))
    nLL <- -sum(dnorm(log(R/Rpred), -0.5 * sigmaR^2, sigmaR, log = TRUE))
    return(nLL)
  }

  if(fix_h) {
    opt <- optimize(solve_SR_par, interval = c(-10, 10), h = h)$minimum
    R0 <- exp(opt)
  } else {
    opt <- nlminb(c(10, 10), solve_SR_par)
    R0 <- exp(opt$par[1])
    if(SR == "BH") h <- 0.2 + 0.8 * ilogit(opt$par[2])
    if(SR == "Ricker") h <- 0.2 + exp(opt$par[2])
  }

  # Virgin reference points
  N0 <- R0 * sum(NPR0)
  E0 <- R0 * EPR0
  VB0 <- R0 * sum(NPR0 * weight * vul)
  B0 <- R0 * sum(NPR0 * weight)

  # Steepness
  if(SR == "BH") {
    Arec <- 4*h/(1-h)/EPR0
    Brec <- (5*h-1)/(1-h)/E0
  }
  if(SR == "Ricker") {
    Arec <- 1/EPR0 * (5*h)^1.25
    Brec <- 1.25 * log(5*h) / E0
  }

  refpt_unfished <- list(h = h, Arec = Arec, Brec = Brec, E0 = E0, R0 = R0, N0 = N0, VB0 = VB0, B0 = B0, EPR0 = EPR0, NPR0 = NPR0)
  refpt_MSY <- SCA_MSY_calc(Arec, Brec, M, weight, mat, vul, SR = SR)
  return(c(refpt_unfished, refpt_MSY))
}
