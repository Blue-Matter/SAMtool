#' @rdname SCA
#' @useDynLib MSEtool
#' @export
SCA2 <- function(x = 1, Data, SR = c("BH", "Ricker"), vulnerability = c("logistic", "dome"), CAA_dist = c("multinomial", "lognormal"),
                 CAA_multiplier = 50, I_type = c("B", "VB", "SSB"), rescale = "mean1", max_age = Data@MaxAge,
                 start = NULL, fix_h = TRUE, fix_F_equilibrium = TRUE, fix_omega = TRUE, fix_sigma = FALSE, fix_tau = TRUE,
                 common_dev = "comp50", integrate = FALSE, silent = TRUE, opt_hess = FALSE, n_restart = ifelse(opt_hess, 0, 1),
                 control = list(iter.max = 2e5, eval.max = 4e5), inner.control = list(),  ...) {
  dependencies <- "Data@Cat, Data@Ind, Data@Mort, Data@L50, Data@L95, Data@CAA, Data@vbK, Data@vbLinf, Data@vbt0, Data@wla, Data@wlb, Data@MaxAge"
  dots <- list(...)
  start <- lapply(start, eval, envir = environment())

  max_age <- as.integer(min(max_age, Data@MaxAge))
  vulnerability <- match.arg(vulnerability)
  CAA_dist <- match.arg(CAA_dist)
  SR <- match.arg(SR)
  I_type <- match.arg(I_type)
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
  est_rec_dev <- rep(1L, length(CAA_n_nominal))
  est_early_rec_dev <- rep(1L, max_age - 1)

  if(rescale == "mean1") rescale <- 1/mean(C_hist)
  data <- list(model = "SCA2", C_hist = C_hist, rescale = rescale, I_hist = I_hist, CAA_hist = t(apply(CAA_hist, 1, function(x) x/sum(x))),
               CAA_n = CAA_n_rescale, n_y = n_y, max_age = max_age, M = M, weight = Wa, mat = mat_age,
               vul_type = vulnerability, I_type = I_type, CAA_dist = CAA_dist, est_early_rec_dev = est_early_rec_dev,
               est_rec_dev = est_rec_dev, yindF = as.integer(0.5 * n_y))
  data$CAA_hist[data$CAA_hist < 1e-8] <- 1e-8

  # Starting values
  params <- list()
  if(!is.null(start)) {
    if(!is.null(start$meanR) && is.numeric(start$meanR)) params$meanRx <- log(start$meanR)
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

  if(is.null(params$meanRx)) {
    params$meanRx <- ifelse(is.null(Data@OM$R0[x]), log(mean(data$C_hist)) + 2, log(1.5 * rescale * Data@OM$R0[x]))
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
    sigmaI <- max(0.05, sdconv(1, Data@CV_Ind[x]), na.rm = TRUE)
    params$log_sigma <- log(sigmaI)
  }
  if(is.null(params$log_tau)) params$log_tau <- log(1)
  params$log_early_rec_dev <- rep(0, max_age - 1)
  params$log_rec_dev <- rep(0, n_y)

  info <- list(Year = Year, data = data, params = params, LH = LH, SR = SR, control = control,
               inner.control = inner.control, fix_h = fix_h)

  map <- list()
  if(any(info$data$C_hist <= 0)) {
    ind <- info$data$C_hist <= 0
    info$params$logF[ind] <- -20
    map_logF <- length(params$logF)
    map_logF[ind] <- NA
    map_logF[!ind] <- 1:sum(!ind)
    map$logF <- factor(map_logF)
  }
  if(fix_F_equilibrium) map$F_equilibrium <- factor(NA)
  if(fix_omega) map$log_omega <- factor(NA)
  if(fix_sigma) map$log_sigma <- factor(NA)
  if(fix_tau) map$log_tau <- factor(NA)
  if(vulnerability == "dome") map$vul_par <- factor(c(1, 2, NA, 3))

  if(is.character(common_dev) && common_dev == "comp50") {
    CAA_all <- colSums(CAA_hist, na.rm = TRUE)/max(colSums(CAA_hist, na.rm = TRUE)) #CAA with max = 1
    CAA_mode <- which.max(CAA_all)[1] # Find max
    comp50_ind <- which(CAA_all[1:CAA_mode] <= 0.5)[1]
    comp50_ind <- comp50_ind[length(comp50_ind)]
    common_dev <- ifelse(is.na(comp50_ind), 0, comp50_ind)
  }
  map_log_rec_dev <- 1:length(params$log_rec_dev)
  if(is.numeric(common_dev) && !is.na(common_dev) && common_dev > 0) {
    ind <- (length(map_log_rec_dev) - common_dev + 1):length(map_log_rec_dev)
    map_log_rec_dev[ind] <- map_log_rec_dev[ind[1]-1]
    map$log_rec_dev <- factor(map_log_rec_dev)
  }

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

  Dev <- c(rev(report$log_early_rec_dev), report$log_rec_dev)

  nll_report <- ifelse(is.character(opt), ifelse(integrate, NA, report$nll), opt$objective)
  Assessment <- new("Assessment", Model = "SCA2", Name = Data@Name, conv = !is.character(SD) && SD$pdHess,
                    FMort = structure(report$F, names = Year),
                    B = structure(report$B, names = Yearplusone),
                    SSB = structure(report$E, names = Yearplusone),
                    VB = structure(report$VB, names = Yearplusone),
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
                    Dev = structure(Dev, names = YearDev),
                    Dev_type = "log-Recruitment deviations",
                    NLL = structure(c(nll_report, report$nll_comp, report$penalty + report$prior),
                                    names = c("Total", "Index", "CAA", "Catch", "Dev", "Penalty")),
                    info = info, obj = obj, opt = opt, SD = SD, TMB_report = report,
                    dependencies = dependencies)


  if(Assessment@conv) {
    info$h <- ifelse(fix_h, Data@steep[x], NA)
    refpt <- SCA_refpt_calc(E = report$E[1:(length(report$E) - 1)], R = report$R[2:length(report$R)], weight = Wa,
                            mat = mat_age, M = M, vul = report$vul, SR = SR, fix_h = fix_h, h = info$h)

    report <- c(report, refpt)
    if(integrate) {
      SE_Early <- sqrt(SD$diag.cov.random[names(SD$par.random) == "log_early_rec_dev"])
      SE_Main <- sqrt(SD$diag.cov.random[names(SD$par.random) == "log_rec_dev"])[map_log_rec_dev]
    } else {
      SE_Early <- sqrt(diag(SD$cov.fixed)[names(SD$par.fixed) == "log_early_rec_dev"])
      SE_Main <- sqrt(diag(SD$cov.fixed)[names(SD$par.fixed) == "log_rec_dev"])[map_log_rec_dev]
    }
    SE_Dev <- c(rev(SE_Early), SE_Main)

    Assessment@FMSY <- report$FMSY
    Assessment@MSY <- report$MSY
    Assessment@BMSY <- report$BMSY
    Assessment@SSBMSY <- report$EMSY
    Assessment@VBMSY <- report$VBMSY
    Assessment@B0 <- report$B0
    Assessment@R0 <- report$R0
    Assessment@N0 <- report$N0
    Assessment@SSB0 <- report$E0
    Assessment@VB0 <- report$VB0
    Assessment@h <- report$h
    Assessment@F_FMSY <- structure(report$F/report$FMSY, names = Year)
    Assessment@B_BMSY <- structure(report$B/report$BMSY, names = Yearplusone)
    Assessment@B_B0 <- structure(report$B/report$B0, names = Yearplusone)
    Assessment@SSB_SSBMSY <- structure(report$E/report$EMSY, names = Yearplusone)
    Assessment@SSB_SSB0 <- structure(report$E/report$E0, names = Yearplusone)
    Assessment@VB_VBMSY <- structure(report$VB/report$VBMSY, names = Yearplusone)
    Assessment@VB_VB0 <- structure(report$VB/report$VB0, names = Yearplusone)
    Assessment@SE_Dev <- structure(SE_Dev, names = YearDev)
    Assessment@TMB_report <- report
  }
  return(Assessment)
}
class(SCA2) <- "Assess"

