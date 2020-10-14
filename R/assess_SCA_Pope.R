
#' @rdname SCA
#' @useDynLib MSEtool
#' @export
SCA_Pope <- function(x = 1, Data, SR = c("BH", "Ricker"), vulnerability = c("logistic", "dome"), CAA_dist = c("multinomial", "lognormal"),
                     CAA_multiplier = 50, I_type = c("B", "VB", "SSB"), rescale = "mean1", max_age = Data@MaxAge,
                     start = NULL, fix_h = TRUE, fix_U_equilibrium = TRUE, fix_sigma = FALSE, fix_tau = TRUE,
                     early_dev = c("comp_onegen", "comp", "all"), late_dev = "comp50", integrate = FALSE,
                     silent = TRUE, opt_hess = FALSE, n_restart = ifelse(opt_hess, 0, 1),
                     control = list(iter.max = 2e5, eval.max = 4e5), inner.control = list(), ...) {
  dependencies <- "Data@Cat, Data@Ind, Data@Mort, Data@L50, Data@L95, Data@CAA, Data@vbK, Data@vbLinf, Data@vbt0, Data@wla, Data@wlb, Data@MaxAge"
  dots <- list(...)
  start <- lapply(start, eval, envir = environment())

  max_age <- as.integer(min(c(max_age, Data@MaxAge)))
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
  data <- list(model = "SCA_Pope", C_hist = C_hist, rescale = rescale, I_hist = I_hist,
               CAA_hist = t(apply(CAA_hist, 1, function(x) x/sum(x))),
               CAA_n = CAA_n_rescale, n_y = n_y, max_age = max_age, M = M,
               weight = Wa, mat = mat_age, vul_type = vulnerability, I_type = I_type,
               SR_type = SR, CAA_dist = CAA_dist, est_early_rec_dev = est_early_rec_dev, est_rec_dev = est_rec_dev)
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
    if(!is.null(start$U_equilibrium) && is.numeric(start$U_equilibrium)) params$U_equilibrium <- start$U_equilibrium
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
  if(is.null(params$U_equilibrium)) params$U_equilibrium <- 0
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
  if(is.null(params$log_sigma)) {
    sigmaI <- max(0.05, sdconv(1, Data@CV_Ind[x]), na.rm = TRUE)
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
  if(fix_h) map$transformed_h <- factor(NA)
  if(fix_U_equilibrium) map$U_equilibrium <- factor(NA)
  if(fix_sigma) map$log_sigma <- factor(NA)
  if(fix_tau) map$log_tau <- factor(NA)
  if(any(!est_early_rec_dev)) map$log_early_rec_dev <- factor(ifelse(est_early_rec_dev, 1:sum(est_early_rec_dev), NA))
  if(any(!est_rec_dev)) map$log_rec_dev <- factor(ifelse(est_rec_dev, 1:sum(est_rec_dev), NA))
  if(vulnerability == "dome") map$vul_par <- factor(c(1, 2, NA, 3))

  random <- NULL
  if(integrate) random <- c("log_early_rec_dev", "log_rec_dev")

  obj <- MakeADFun(data = info$data, parameters = info$params, hessian = TRUE,
                   map = map, random = random, DLL = "MSEtool", inner.control = inner.control, silent = silent)

  # Add starting values for rec-devs and increase R0 start value if U is too high (> 0.975)
  high_U <- try(obj$report(c(obj$par, obj$env$last.par[obj$env$random]))$penalty > 0, silent = TRUE)
  if(!is.character(high_U) && high_U) {
    Recruit <- try(Data@Rec[x, ], silent = TRUE)
    if(is.numeric(Recruit) && length(Recruit) == n_y && any(!is.na(Recruit))) {
      log_rec_dev <- log(Recruit/mean(Recruit, na.rm = TRUE))
      log_rec_dev[is.na(est_rec_dev) | is.na(log_rec_dev) | is.infinite(log_rec_dev)] <- 0
      info$params$log_rec_dev <- log_rec_dev

      obj <- MakeADFun(data = info$data, parameters = info$params, hessian = TRUE,
                       map = map, random = random, DLL = "MSEtool", inner.control = inner.control, silent = silent)
    }
    while(obj$par["R0x"] < 30 && obj$report(c(obj$par, obj$env$last.par[obj$env$random]))$penalty > 0) {
      obj$par["R0x"] <- obj$par["R0x"] + 1
    }
  }

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
  Assessment <- new("Assessment", Model = "SCA_Pope", Name = Data@Name, conv = !is.character(SD) && SD$pdHess,
                    B0 = report$B0, R0 = report$R0, N0 = report$N0,
                    SSB0 = report$E0, VB0 = report$VB0,
                    h = report$h, U = structure(report$U, names = Year),
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
                    Catch = structure(colSums(t(report$CAApred) * Wa), names = Year),
                    Index = structure(report$Ipred, names = Year),
                    C_at_age = report$CAApred,
                    Dev = Dev, Dev_type = "log-Recruitment deviations",
                    NLL = structure(c(nll_report, report$nll_comp, report$penalty),
                                    names = c("Total", "Index", "CAA", "Dev", "Penalty")),
                    info = info, obj = obj, opt = opt, SD = SD, TMB_report = report,
                    dependencies = dependencies)

  if(Assessment@conv) {
    ref_pt <- SCA_Pope_MSY_calc(Arec = report$Arec, Brec = report$Brec, M = M, weight = Wa, mat = mat_age, vul = report$vul, SR = SR)
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

    Assessment@UMSY <- report$UMSY
    Assessment@MSY <- report$MSY
    Assessment@BMSY <- report$BMSY
    Assessment@SSBMSY <- report$EMSY
    Assessment@VBMSY <- report$VBMSY
    Assessment@U_UMSY <- structure(report$U/report$UMSY, names = Year)
    Assessment@B_BMSY <- structure(report$B/report$BMSY, names = Yearplusone)
    Assessment@SSB_SSBMSY <- structure(report$E/report$EMSY, names = Yearplusone)
    Assessment@VB_VBMSY <- structure(report$VB/report$VBMSY, names = Yearplusone)
    Assessment@Dev <- Dev
    Assessment@SE_Dev <- SE_Dev
    Assessment@TMB_report <- report
  }
  return(Assessment)
}
class(SCA_Pope) <- "Assess"



SCA_Pope_MSY_calc <- function(Arec, Brec, M, weight, mat, vul, SR = c("BH", "Ricker")) {
  SR <- match.arg(SR)
  maxage <- length(M)

  solveMSY <- function(logit_U) {
    U <- ilogit(logit_U)
    surv <- exp(-M) * (1 - vul * U)
    NPR <- c(1, cumprod(surv[1:(maxage-1)]))
    NPR[maxage] <- NPR[maxage]/(1 - surv[maxage])
    EPR <- sum(NPR * mat * weight)
    if(SR == "BH") Req <- (Arec * EPR - 1)/(Brec * EPR)
    if(SR == "Ricker") Req <- log(Arec * EPR)/(Brec * EPR)
    CPR <- vul * U * NPR * exp(-0.5 * M)
    Yield <- Req * sum(CPR * weight)
    return(-1 * Yield)
  }

  opt2 <- optimize(solveMSY, interval = c(logit(0.01), logit(0.99)))
  UMSY <- ilogit(opt2$minimum)
  MSY <- -1 * opt2$objective
  VBMSY <- MSY/UMSY

  surv_UMSY <- exp(-M) * (1 - vul * UMSY)
  NPR_UMSY <- c(1, cumprod(surv_UMSY[1:(maxage-1)]))
  NPR_UMSY[maxage] <- NPR_UMSY[maxage]/(1 - surv_UMSY[maxage])

  RMSY <- VBMSY/sum(vul * NPR_UMSY * weight)
  BMSY <- RMSY * sum(NPR_UMSY * weight)
  EMSY <- RMSY * sum(NPR_UMSY * weight * mat)

  return(list(UMSY = UMSY, MSY = MSY, VBMSY = VBMSY, RMSY = RMSY, BMSY = BMSY, EMSY = EMSY))
}


vul_fn <- function(vul_par, maxage, type) {
  age <- 1:maxage

  if(type == "logistic") {
    a50 <- vul_par[1]
    a95 <- a50 + exp(vul_par[2])
    vul <- 1/(1 + exp(-log(19) * (age - a50)/(a95 - a50)))
  }
  if(type == "dome") {
    sd_asc <- exp(vul_par[1])
    mu_asc <- vul_par[2]
    mu_des <- mu_asc + exp(vul_par[3])
    sd_des <- exp(vul_par[4])

    denom_asc <- dnorm(mu_asc, mu_asc, sd_asc)
    denom_des <- dnorm(mu_des, mu_des, sd_des)

    vul <- rep(NA, maxage)
    for(i in age) {
      if(i <= mu_asc) {
        vul[i] <- dnorm(i, mu_asc, sd_asc)/denom_asc
      } else if(i <= mu_des) {
        vul[i] <- 1
      } else {
        vul[i] <- dnorm(i, mu_des, sd_des)/denom_des
      }
    }
  }
  return(vul)
}

