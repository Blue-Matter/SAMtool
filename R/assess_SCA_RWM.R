
#' SCA with random walk in M
#'
#' \code{SCA_RWM} is a modification of \link{SCA} that incorporates a random walk in M (constant with age. 
#' Set the variance to a small value (0.01) in order to estimate a time and age-constant M.
#' 
#' @param x A position in the Data object (by default, equal to one for assessments).
#' @param Data An object of class Data
#' @param SR Stock-recruit function (either \code{"BH"} for Beverton-Holt or \code{"Ricker"}).
#' @param vulnerability Whether estimated vulnerability is \code{"logistic"} or \code{"dome"} (double-normal).
#' See details for parameterization.
#' @param CAA_dist Whether a multinomial or lognormal distribution is used for likelihood of the catch-at-age matrix. See details.
#' @param CAA_multiplier Numeric for data weighting of catch-at-age matrix if \code{CAA_hist = "multinomial"}. Otherwise ignored. See details.
#' @param I_type Whether the index surveys population biomass (B; this is the default in the MSEtool operating model),
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
#' The model estimates year-specific M (constant with age) as a random walk.
#' 
#' The starting value for the first year M (M_start) is \code{Data@@Mort[x]}, and the fixed SD of the random walk
#' in logspace (\code{tau_M}) is 0.05. Alternative values can be provided in the start list (see example). 
#' 
#' See \link{SCA} for all other information about the structure and setup of the model.
#' 
#' The SCA builds in a stock-recruit relationship into the model. Annual unfished andMSY reference points are still 
#' calculated and reported in TMB_report.
#'
#' @examples
#' res <- SCA_RWM(Data = MSEtool::SimulatedData, start = list(M_start = 0.4, tau_M = 0.05))
#' res2 <- SCA(Data = MSEtool::SimulatedData)
#' res3 <- SCA_RWM(Data = MSEtool::SimulatedData, start = list(M_start = 0.4, tau_M = 0.001))
#'
#' compare_models(res, res2, res3)
#'
#' \dontrun{
#' plot(res)
#' }
#' @author Q. Huynh
#' @return An object of class \linkS4class{Assessment}.
#' @seealso \link{SCA}
#' @export
SCA_RWM <- function(x = 1, Data, SR = c("BH", "Ricker"), vulnerability = c("logistic", "dome"), CAA_dist = c("multinomial", "lognormal"),
                    CAA_multiplier = 50, I_type = c("B", "VB", "SSB"), rescale = "mean1", max_age = Data@MaxAge,
                    start = NULL, fix_h = TRUE, fix_F_equilibrium = TRUE, fix_omega = TRUE, fix_sigma = FALSE, fix_tau = TRUE,
                    early_dev = c("comp_onegen", "comp", "all"), late_dev = "comp50", integrate = FALSE,
                    silent = TRUE, opt_hess = FALSE, n_restart = ifelse(opt_hess, 0, 1),
                    control = list(iter.max = 2e5, eval.max = 4e5), inner.control = list(), ...) {
  dependencies <- "Data@Cat, Data@Ind, Data@Mort, Data@L50, Data@L95, Data@CAA, Data@vbK, Data@vbLinf, Data@vbt0, Data@wla, Data@wlb, Data@MaxAge"
  dots <- list(...)
  start <- lapply(start, eval, envir = environment())

  max_age <- as.integer(min(max_age, Data@MaxAge))
  n_age <- max_age + 1
  vulnerability <- match.arg(vulnerability)
  CAA_dist <- match.arg(CAA_dist)
  SR <- match.arg(SR)
  I_type <- match.arg(I_type)
  if(is.character(early_dev)) early_dev <- match.arg(early_dev)
  if(is.numeric(early_dev)) stopifnot(early_dev < length(Data@Year))

  if(any(names(dots) == "yind")) {
    yind <- eval(dots$yind)
  } else {
    yind <- which(!is.na(Data@Cat[x, ]))[1]
    yind <- yind:length(Data@Cat[x, ])
  }
  Year <- Data@Year[yind]
  C_hist <- Data@Cat[x, yind]
  if(any(is.na(C_hist) | C_hist < 0)) warning("Error. Catch time series is not complete.")
  if(I_type == "B") {
    I_hist <- Data@Ind[x, yind]
  } else if(I_type == "VB") {
    I_hist <- Data@VInd[x, yind]
  } else {
    I_hist <- Data@Sp_Ind[x, yind]
  }
  Data <- expand_comp_matrix(Data, "CAA") # Make sure dimensions of CAA match that in catch (nyears).
  CAA_hist <- Data@CAA[x, yind, 1:n_age]
  if(max_age < Data@MaxAge) CAA_hist[, n_age] <- rowSums(Data@CAA[x, yind, n_age:(Data@MaxAge+1)], na.rm = TRUE)

  CAA_n_nominal <- rowSums(CAA_hist)
  if(CAA_multiplier <= 1) {
    CAA_n_rescale <- CAA_multiplier * CAA_n_nominal
  } else CAA_n_rescale <- pmin(CAA_multiplier, CAA_n_nominal)

  n_y <- length(C_hist)
  #M <- rep(Data@Mort[x], n_age)
  a <- Data@wla[x]
  b <- Data@wlb[x]
  Linf <- Data@vbLinf[x]
  K <- Data@vbK[x]
  t0 <- Data@vbt0[x]
  La <- Linf * (1 - exp(-K * (c(0:max_age) - t0)))
  Wa <- a * La ^ b
  A50 <- min(0.5 * max_age, iVB(t0, K, Linf, Data@L50[x]))
  A95 <- max(A50+0.5, iVB(t0, K, Linf, Data@L95[x]))
  mat_age <- c(0, 1/(1 + exp(-log(19) * (c(1:max_age) - A50)/(A95 - A50)))) # Age-0 is immature
  mat_age <- mat_age/max(mat_age)
  LH <- list(LAA = La, WAA = Wa, Linf = Linf, K = K, t0 = t0, a = a, b = b, A50 = A50, A95 = A95)

  if(early_dev == "all") {
    est_early_rec_dev <- rep(1, n_age - 1)
    est_rec_dev <- rep(1, n_y)
  } else if(early_dev == "comp") {
    est_early_rec_dev <- rep(0, n_age-1)
    ind1 <- which(!is.na(CAA_n_nominal))[1]
    est_rec_dev <- ifelse(1:n_y < ind1, 0, 1)
  } else if(early_dev == "comp_onegen") {
    ind1 <- which(!is.na(CAA_n_nominal))[1] - n_age
    if(ind1 < 0) {
      early_start <- n_age + ind1
      est_early_rec_dev <- rev(ifelse(c(1:(n_age-1)) < early_start, 0, 1))
      est_rec_dev <- rep(1, n_y)
    } else {
      est_early_rec_dev <- rep(0, n_age - 1)
      est_rec_dev <- ifelse(1:n_y < ind1, 0, 1)
    }
  } else if(is.numeric(early_dev)) {
    if(early_dev > 1) {
      est_early_rec_dev <- rep(0, n_age-1)
      est_rec_dev <- ifelse(1:n_y >= early_dev, 1, 0)
    } else {
      ind1 <- early_dev - 1
      est_early_rec_dev <- c(rep(1, ind1), rep(NA, n_age-ind1-1))
      est_rec_dev <- rep(1, n_y)
    }
  } else {
    stop("Invalid early_dev argument.")
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
  data <- list(model = "SCA_RWM", C_hist = C_hist, rescale = rescale, I_hist = I_hist,
               CAA_hist = t(apply(CAA_hist, 1, function(x) x/sum(x))),
               CAA_n = CAA_n_rescale, n_y = n_y, n_age = n_age, #M = M,
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
      params$log_F_dev <- Fstart
    }
    if(!is.null(start$M_start) && is.numeric(start$M_start)) params$log_M_start <- log(start$M_start)
    
    if(!is.null(start$omega) && is.numeric(start$omega)) params$log_omega <- log(start$omega)
    if(!is.null(start$sigma) && is.numeric(start$sigma)) params$log_sigma <- log(start$sigma)
    if(!is.null(start[["tau"]]) && is.numeric(start[["tau"]])) params$log_tau <- log(start[["tau"]])
    if(!is.null(start$tau_M) && is.numeric(start$tau_M)) params$log_tau_M <- log(start$tau_M)
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
  if(is.null(params$log_F_dev)) {
    Fstart <- numeric(n_y)
    Fstart[data$yindF + 1] <- log(0.75 * mean(Data@Mort[x]))
    params$log_F_dev <- Fstart
  }
  if(is.null(params$log_M_start)) params$log_M_start <- log(Data@Mort[x])
  params$log_M_walk <- rep(0, n_y - 1)
  
  if(is.null(params$log_omega)) {
    sigmaC <- max(0.01, sdconv(1, Data@CV_Cat[x]), na.rm = TRUE)
    params$log_omega <- log(sigmaC)
  }
  if(is.null(params$log_sigma)) {
    sigmaI <- max(0.01, sdconv(1, Data@CV_Ind[x]), na.rm = TRUE)
    params$log_sigma <- log(sigmaI)
  }
  if(is.null(params[["log_tau"]])) {
    tau_start <- ifelse(is.na(Data@sigmaR[x]), 0.6, Data@sigmaR[x])
    params$log_tau <- log(tau_start)
  }
  if(is.null(params$log_tau_M)) params$log_tau_M <- log(0.05)
  params$log_M_walk <- rep(0, n_y - 1)
  
  params$log_early_rec_dev <- rep(0, n_age - 1)
  params$log_rec_dev <- rep(0, n_y)

  info <- list(Year = Year, data = data, params = params, LH = LH, control = control,
               inner.control = inner.control)

  map <- list()
  if(any(info$data$C_hist <= 0)) {
    ind <- info$data$C_hist <= 0
    info$params$log_F_dev[ind] <- -20
    map_logF <- length(params$log_F_dev)
    map_logF[ind] <- NA
    map_logF[!ind] <- 1:sum(!ind)
    map$log_F_dev <- factor(map_logF)
  }
  if(fix_h) map$transformed_h <- factor(NA)
  if(fix_F_equilibrium) map$F_equilibrium <- factor(NA)
  if(fix_omega) map$log_omega <- factor(NA)
  if(fix_sigma) map$log_sigma <- factor(NA)
  if(fix_tau) map$log_tau <- factor(NA)
  if(any(!est_early_rec_dev)) map$log_early_rec_dev <- factor(ifelse(est_early_rec_dev, 1:sum(est_early_rec_dev), NA))
  if(any(!est_rec_dev)) map$log_rec_dev <- factor(ifelse(est_rec_dev, 1:sum(est_rec_dev), NA))
  if(vulnerability == "dome") map$vul_par <- factor(c(1, 2, NA, 3))
  map$log_tau_M <- factor(NA)

  random <- NULL
  if(integrate) random <- c("log_early_rec_dev", "log_rec_dev", "log_M_walk")

  obj <- MakeADFun(data = info$data, parameters = info$params, hessian = TRUE,
                   map = map, random = random, DLL = "SAMtool", inner.control = inner.control, silent = silent)

  mod <- optimize_TMB_model(obj, control, opt_hess, n_restart)
  opt <- mod[[1]]
  SD <- mod[[2]]
  report <- obj$report(obj$env$last.par.best)
  report$NPR0 <- do.call(cbind, report$NPR0)
  report$E0_tv <- report$EPR0 * report$R0
  report$E0_tv <- c(report$E0_tv, report$E0_tv[length(report$E0_tv)])

  Yearplusone <- c(Year, max(Year) + 1)
  YearEarly <- (Year[1] - n_age + 1):(Year[1] - 1)
  YearDev <- c(YearEarly, Year)
  YearR <- c(YearDev, max(YearDev) + 1)
  R <- c(rev(report$R_early), report$R)

  Dev <- structure(c(rev(report$log_early_rec_dev), report$log_rec_dev), names = YearDev)
  report$dynamic_SSB0 <- SCA_dynamic_SSB0(obj) %>% structure(names = Yearplusone)
  
  nll_report <- ifelse(is.character(opt), ifelse(integrate, NA, report$nll), opt$objective)
  Assessment <- new("Assessment", Model = "SCA_RWM", Name = Data@Name, conv = !is.character(SD) && SD$pdHess,
                    B0 = report$B0, R0 = report$R0, N0 = report$N0,
                    SSB0 = report$E0, VB0 = report$VB0,
                    h = report$h, FMort = structure(report$F, names = Year),
                    B = structure(report$B, names = Yearplusone),
                    B_B0 = structure(report$B/report$B0, names = Yearplusone),
                    SSB = structure(report$E, names = Yearplusone),
                    SSB_SSB0 = structure(report$E/report$E0_tv, names = Yearplusone),
                    VB = structure(report$VB, names = Yearplusone),
                    VB_VB0 = structure(report$VB/report$VB0, names = Yearplusone),
                    R = structure(R, names = YearR),
                    N = structure(rowSums(report$N), names = Yearplusone),
                    N_at_age = report$N,
                    Selectivity = matrix(report$vul, nrow = length(Year),
                                         ncol = n_age, byrow = TRUE),
                    Obs_Catch = structure(C_hist, names = Year),
                    Obs_Index = structure(I_hist, names = Year),
                    Obs_C_at_age = CAA_hist,
                    Catch = structure(report$Cpred, names = Year),
                    Index = structure(report$Ipred, names = Year),
                    C_at_age = report$CAApred,
                    Dev = Dev, Dev_type = "log-Recruitment deviations",
                    NLL = structure(c(nll_report, report$nll_comp, report$penalty + report$prior),
                                    names = c("Total", "Index", "CAA", "Catch", "RecDev", "M_Walk", "Penalty")),
                    info = info, obj = obj, opt = opt, SD = SD, TMB_report = report,
                    dependencies = dependencies)

  if(Assessment@conv) {
    year_specific_ref_pt <- lapply(report$M, function(x) {
      SCA_MSY_calc(Arec = report$Arec, Brec = report$Brec, M = rep(x, n_age), 
                   weight = Wa, mat = mat_age, vul = report$vul, SR = SR)
    })
    ref_pt <- lapply(names(year_specific_ref_pt[[1]]), function(x) sapply(year_specific_ref_pt, getElement, x)) %>%
      structure(names = names(year_specific_ref_pt[[1]]))
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
    #Assessment@B_BMSY <- structure(report$B/report$BMSY, names = Yearplusone)
    #Assessment@SSB_SSBMSY <- structure(report$E/report$EMSY, names = Yearplusone)
    #Assessment@VB_VBMSY <- structure(report$VB/report$VBMSY, names = Yearplusone)
    Assessment@Dev <- Dev
    Assessment@SE_Dev <- SE_Dev
    Assessment@TMB_report <- report
  }
  return(Assessment)
}
class(SCA_RWM) <- "Assess"
