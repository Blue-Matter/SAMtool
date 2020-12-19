#' @rdname SCA
#' @useDynLib SAMtool
#' @export
SCA2 <- function(x = 1, Data, AddInd = "B", SR = c("BH", "Ricker"), vulnerability = c("logistic", "dome"), CAA_dist = c("multinomial", "lognormal"),
                 CAA_multiplier = 50, rescale = "mean1", max_age = Data@MaxAge, start = NULL, 
                 fix_h = TRUE, fix_F_equilibrium = TRUE, fix_omega = TRUE, fix_tau = TRUE, LWT = NULL,
                 common_dev = "comp50", integrate = FALSE, silent = TRUE, opt_hess = FALSE, n_restart = ifelse(opt_hess, 0, 1),
                 control = list(iter.max = 2e5, eval.max = 4e5), inner.control = list(),  ...) {
  dependencies <- "Data@Cat, Data@Ind, Data@Mort, Data@L50, Data@L95, Data@CAA, Data@vbK, Data@vbLinf, Data@vbt0, Data@wla, Data@wlb, Data@MaxAge"
  dots <- list(...)
  start <- lapply(start, eval, envir = environment())

  max_age <- as.integer(min(max_age, Data@MaxAge))
  n_age <- max_age + 1
  vulnerability <- match.arg(vulnerability)
  CAA_dist <- match.arg(CAA_dist)
  SR <- match.arg(SR)
  
  if(any(names(dots) == "yind")) {
    yind <- eval(dots$yind)
  } else {
    yind <- which(!is.na(Data@Cat[x, ]))[1]
    yind <- yind:length(Data@Cat[x, ])
  }
  Year <- Data@Year[yind]
  C_hist <- Data@Cat[x, yind]
  if(any(is.na(C_hist) | C_hist < 0)) warning("Error. Catch time series is not complete.")
  
  Data <- expand_comp_matrix(Data, "CAA") # Make sure dimensions of CAA match that in catch (nyears).
  CAA_hist <- Data@CAA[x, yind, 1:n_age]
  if(max_age < Data@MaxAge) CAA_hist[, n_age] <- rowSums(Data@CAA[x, yind, n_age:(Data@MaxAge+1)], na.rm = TRUE)

  CAA_n_nominal <- rowSums(CAA_hist)
  if(CAA_multiplier <= 1) {
    CAA_n_rescale <- CAA_multiplier * CAA_n_nominal
  } else CAA_n_rescale <- pmin(CAA_multiplier, CAA_n_nominal)

  n_y <- length(C_hist)
  M <- rep(Data@Mort[x], n_age)
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
  est_rec_dev <- rep(1L, length(CAA_n_nominal))
  est_early_rec_dev <- rep(1L, n_age - 1)

  if(rescale == "mean1") rescale <- 1/mean(C_hist)
  
  Ind <- lapply(AddInd, Assess_I_hist, Data = Data, x = x, yind = yind)
  I_hist <- do.call(cbind, lapply(Ind, getElement, "I_hist"))
  I_sd <- do.call(cbind, lapply(Ind, getElement, "I_sd")) %>% pmax(0.05)
  I_units <- do.call(c, lapply(Ind, getElement, "I_units"))
  I_vul <- vapply(AddInd, function(xx) {
    if(xx == "B") {
      return(rep(1, n_age))
    } else if(xx == "SSB") {
      return(mat_age)
    } else if(xx == "VB") {
      return(rep(0, n_age))
    } else {
      return(Data@AddIndV[x, suppressWarnings(as.numeric(xx)), 1:n_age])
    }
  }, numeric(n_age))
  nsurvey <- ncol(I_hist)
  if(is.null(LWT)) LWT <- rep(1, nsurvey)
  if(length(LWT) != nsurvey) stop("LWT needs to be a vector of length ", nsurvey)
  
  data <- list(model = "SCA2", C_hist = C_hist, rescale = rescale, I_hist = I_hist, 
               I_sd = I_sd, I_units = I_units, I_vul = I_vul, abs_I = rep(0, nsurvey), nsurvey = nsurvey, LWT = LWT,
               CAA_hist = t(apply(CAA_hist, 1, function(x) x/sum(x))),
               CAA_n = CAA_n_rescale, n_y = n_y, n_age = n_age, M = M, weight = Wa, mat = mat_age,
               vul_type = vulnerability, CAA_dist = CAA_dist, est_early_rec_dev = est_early_rec_dev,
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
      params$log_F_dev <- Fstart
    }

    if(!is.null(start$omega) && is.numeric(start$omega)) params$log_omega <- log(start$omega)
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
  if(is.null(params$log_F_dev)) {
    Fstart <- numeric(n_y)
    Fstart[data$yindF + 1] <- log(0.75 * mean(data$M))
    params$log_F_dev <- Fstart
  }
  if(is.null(params$log_omega)) {
    sigmaC <- max(0.01, sdconv(1, Data@CV_Cat[x]), na.rm = TRUE)
    params$log_omega <- log(sigmaC)
  }
  if(is.null(params$log_tau)) params$log_tau <- log(1)
  params$log_early_rec_dev <- rep(0, n_age - 1)
  params$log_rec_dev <- rep(0, n_y)

  info <- list(Year = Year, data = data, params = params, LH = LH, SR = SR, control = control,
               inner.control = inner.control, fix_h = fix_h)

  map <- list()
  if(any(info$data$C_hist <= 0)) {
    ind <- info$data$C_hist <= 0
    info$params$log_F_dev[ind] <- -20
    map_logF <- length(params$log_F_dev)
    map_logF[ind] <- NA
    map_logF[!ind] <- 1:sum(!ind)
    map$log_F_dev <- factor(map_logF)
  }
  if(fix_F_equilibrium) map$F_equilibrium <- factor(NA)
  if(fix_omega) map$log_omega <- factor(NA)
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
                   map = map, random = random, DLL = "SAMtool", inner.control = inner.control, silent = silent)

  mod <- optimize_TMB_model(obj, control, opt_hess, n_restart)
  opt <- mod[[1]]
  SD <- mod[[2]]

  report <- obj$report(obj$env$last.par.best)

  Yearplusone <- c(Year, max(Year) + 1)
  YearEarly <- (Year[1] - n_age + 1):(Year[1] - 1)
  YearDev <- c(YearEarly, Year)
  YearR <- c(YearDev, max(YearDev) + 1)
  R <- c(rev(report$R_early), report$R)

  Dev <- c(rev(report$log_early_rec_dev), report$log_rec_dev)
  report$dynamic_SSB0 <- SCA_dynamic_SSB0(obj) %>% structure(names = Yearplusone)

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
                                         ncol = n_age, byrow = TRUE),
                    Obs_Catch = structure(C_hist, names = Year),
                    Obs_Index = structure(I_hist, dimnames = list(Year, paste0("Index_", 1:nsurvey))),
                    Obs_C_at_age = CAA_hist,
                    Catch = structure(report$Cpred, names = Year),
                    Index = structure(report$Ipred, dimnames = list(Year, paste0("Index_", 1:nsurvey))),
                    C_at_age = report$CAApred,
                    Dev = structure(Dev, names = YearDev),
                    Dev_type = "log-Recruitment deviations",
                    NLL = structure(c(nll_report, report$nll_comp, report$prior, report$penalty),
                                    names = c("Total", paste0("Index_", 1:nsurvey), "CAA", "Catch", "Dev", "Prior", "Penalty")),
                    info = info, obj = obj, opt = opt, SD = SD, TMB_report = report,
                    dependencies = dependencies)

  if(Assessment@conv) {
    info$h <- ifelse(fix_h, Data@steep[x], NA)
    ref_pt <- ref_pt_SCA2(E = report$E, R = report$R, weight = Wa,
                          mat = mat_age, M = M, vul = report$vul, SR = SR, fix_h = fix_h, h = info$h)

    report <- c(report, ref_pt[-length(ref_pt)])
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
    
    catch_eq <- function(Ftarget) {
      projection_SCA2(Assessment, Ftarget = Ftarget, p_years = 1, p_sim = 1, obs_error = list(matrix(1, 1, 1), matrix(1, 1, 1)),
                      process_error = matrix(1, 1, 1)) %>% slot("Catch") %>% as.vector()
    }
    Assessment@forecast <- list(per_recruit = ref_pt[[length(ref_pt)]], catch_eq = catch_eq)
  }
  return(Assessment)
}
class(SCA2) <- "Assess"


ref_pt_SCA2 <- function(E, R, weight, mat, M, vul, SR, fix_h, h) {
  # Unfished spawners per recruit
  n_age <- length(M)
  surv0 <- exp(-M)
  NPR0 <- c(1, cumprod(surv0[1:(n_age-1)]))
  NPR0[n_age] <- NPR0[n_age]/(1 - exp(-M[n_age]))
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
  
  # Unfished reference points
  N0 <- R0 * sum(NPR0)
  E0 <- R0 * EPR0
  VB0 <- R0 * sum(NPR0 * weight * vul)
  B0 <- R0 * sum(NPR0 * weight)
  
  # Alpha/Beta from steepness
  if(SR == "BH") {
    Arec <- 4*h/(1-h)/EPR0
    Brec <- (5*h-1)/(1-h)/E0
  }
  if(SR == "Ricker") {
    Arec <- 1/EPR0 * (5*h)^1.25
    Brec <- 1.25 * log(5*h) / E0
  }
  
  refpt_unfished <- list(h = h, Arec = Arec, Brec = Brec, E0 = E0, R0 = R0, N0 = N0, VB0 = VB0, B0 = B0, EPR0 = EPR0, NPR0 = NPR0)
  refpt_MSY <- ref_pt_SCA(Arec, Brec, M, weight, mat, vul, SR = SR)
  return(c(refpt_unfished, refpt_MSY))
}
