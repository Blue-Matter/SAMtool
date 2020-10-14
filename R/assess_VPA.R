#' Virtual population analysis (VPA)
#'
#' A VPA model that back-calculates abundance-at-age assuming that the catch-at-age is known without error and tuned to an index.
#' The population dynamics equations are primarily drawn from VPA-2BOX (Porch 2018). MSY reference points are then calculated from the
#' VPA output.
#'
#' @param x A position in the Data object (by default, equal to one for assessments).
#' @param Data An object of class Data
#' @param expanded Whether the catch at age in \code{Data} has been expanded. If \code{FALSE}, then the catch in weight
#' should be provided in \code{Data@@Cat} so that the function can calculate annual expansion factors.
#' @param SR Stock-recruit function (either \code{"BH"} for Beverton-Holt or \code{"Ricker"}) for calculating MSY reference points.
#' @param vulnerability Whether the terminal year vulnerability is \code{"logistic"} or \code{"dome"} (double-normal). If \code{"free"},
#' independent F's are calculated in the terminal year (subject to the assumed ratio of F of the plus-group to the previous age class).
#' See details for parameterization.
#' @param I_type Whether the index surveys population biomass (B; this is the default in the DLMtool operating model),
#' vulnerable biomass (VB), or spawning stock biomass (SSB).
#' @param start Optional list of starting values. Entries can be expressions that are evaluated in the function. See details.
#' @param fix_Fratio Logical, whether the ratio of F of the plus-group to the previous age class is fixed in the model.
#' @param fix_h Logical, whether to fix steepness to value in \code{Data@@steep}. This only affects
#' calculation of reference points.
#' @param fix_sigma Logical, whether the standard deviation of the index is fixed. If \code{TRUE},
#' sigma is fixed to value provided in \code{start} (if provided), otherwise, value based on \code{Data@@CV_Ind}.
#' @param vul_pen A length two vector that parameterizes how the model constrains the vulnerability in the most recent years. The first number
#' is the number of years in which vulnerability will be constrained (as a random walk), the second number is the standard deviation of the random walk.
#' @param R_pen A length two vector that parameterizes how the model constrains the recruitment in the most recent years. The first number
#' is the number of years in which recruitment will be constrained (as a random walk), the second number is the standard deviation of the random walk.
#' @param nitF The number of iterations for solving F in the model (via Newton's method).
#' @param silent Logical, passed to \code{\link[TMB]{MakeADFun}}, whether TMB
#' will print trace information during optimization. Used for dignostics for model convergence.
#' @param opt_hess Logical, whether the hessian function will be passed to \code{\link[stats]{nlminb}} during optimization
#' (this generally reduces the number of iterations to convergence, but is memory and time intensive and does not guarantee an increase
#' in convergence rate). Ignored if \code{integrate = TRUE}.
#' @param n_restart The number of restarts (calls to \code{\link[stats]{nlminb}}) in the optimization procedure, so long as the model
#' hasn't converged. The optimization continues from the parameters from the previous (re)start.
#' @param control A named list of agruments for optimization to be passed to
#' \code{\link[stats]{nlminb}}.
#' @param ... Other arguments to be passed.
#' @details
#' The VPA is initialized by estimating the terminal F-at-age. Parameter \code{F_term} is the apical terminal F if
#' a functional form for vulnerability is used in the terminal year. If the terminal F-at-age are otherwise independent parameters,
#' \code{F_term} is the F for the reference age which is half the maximum age. Once terminal-year abundance is
#' estimated, the abundance in historical years can be back-calculated. The oldest age group is a plus-group, and requires
#' an assumption regarding the ratio of F's between the plus-group and the next youngest age class. The F-ratio can
#' be fixed (default) or estimated.
#'
#'
#'
#' For \code{start} (optional), a named list of starting values of estimates can be provided for:
#' \itemize{
#' \item \code{F_term} The terminal year fishing mortality.
#' \item \code{F_ratio} The ratio of F in the plus-group to the next youngest age. If not provided, a value of 1 is used.
#' \item \code{vul_par} Vulnerability parameters in the terminal year. This will be of length 2 vector for \code{"logistic"} or length 4 for
#' \code{"dome"}, see \link{SCA} for further documentation on parameterization. For option \code{"free"}, this will be a vector of length
#' A-2 where A is the number of age classes in the model. To estimate parameters, vulnerability is initially set to one at half the max age
#' (and subsequently re-calculated relative to the maximum F experienced in that year). Vulnerability in the plus-group is also constrained
#' by the Fratio.
#' \item \code{sigma} Standard deviation of the index. If not provided, the value based on \code{Data@@CV_Ind} is used.
#' }
#' @return An object of class \linkS4class{Assessment}. The F vector is the apical fishing mortality experienced by any
#' age class in a given year. The U vector is the ratio of catch (weight) and vulnerable biomass, which may be a better
#' description of fishing pressure (and UMSY = MSY/VBMSY).
#' @references
#' Porch, C.E. 2018. VPA-2BOX 4.01 User Guide. NOAA Tech. Memo. NMFS-SEFSC-726. 67 pp.
#' @export
VPA <- function(x = 1, Data, expanded = FALSE, SR = c("BH", "Ricker"), vulnerability = c("logistic", "dome", "free"),
                I_type = c("B", "VB", "SSB"), start = NULL, fix_h = TRUE,
                fix_sigma = FALSE, fix_Fratio = TRUE, vul_pen = c(3, 0.4), R_pen = c(3, Data@sigmaR[x]), nitF = 5L,
                silent = TRUE, opt_hess = FALSE, n_restart = ifelse(opt_hess, 0, 1),
                control = list(iter.max = 2e5, eval.max = 4e5), ...) {
  dependencies <- "Data@Cat, Data@CAA, Data@Ind, Data@Mort, Data@L50, Data@L95, Data@CAA, Data@vbK, Data@vbLinf, Data@vbt0, Data@wla, Data@wlb, Data@MaxAge"
  dots <- list(...)
  start <- lapply(start, eval, envir = environment())

  vulnerability <- match.arg(vulnerability)
  SR <- match.arg(SR)
  I_type <- match.arg(I_type)

  # Life history
  max_age <- Data@MaxAge
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

  # Data
  Data <- expand_comp_matrix(Data, "CAA") # Make sure dimensions of CAA match that in catch (nyears).

  if(any(names(dots) == "yind")) {
    yind <- eval(dots$yind)
  } else {
    yind <- 1:nrow(Data@CAA[x, , ])
  }

  Year <- Data@Year[yind]
  I_hist <- Data@Ind[x, yind]
  CAA_hist <- Data@CAA[x, yind, ]

  if(!expanded) {
    C_hist <- Data@Cat[x, yind]
    if(any(is.na(C_hist) | C_hist < 0)) warning("Error. Catch time series is not complete.")
    expansion_factors <- C_hist/colSums(t(CAA_hist) * Wa)
    CAA_hist <- CAA_hist * expansion_factors
  }

  if(any(names(dots) == "ages")) {
    ages <- dots$ages

    CAA_hist2 <- CAA_hist[, 1:max(ages)]
    CAA_hist2[, max(ages)] <- rowSums(CAA_hist[, max(ages):max_age, drop = FALSE])
    if(min(ages) > 1) {
      CAA_hist2 <- CAA_hist2[, -c(1:(min(ages)-1)), drop = FALSE]
      CAA_hist2[, 1] <- rowSums(CAA_hist[, 1:min(ages), drop = FALSE])
    }

  } else {

    max_age2 <- max_age  # Reduce max-age until no zero's are observed
    CAA_hist2 <- CAA_hist
    while(any(CAA_hist2[, max_age2] <= 0)) {
      max_age2 <- max_age2 - 1
      CAA_hist2 <- CAA_hist[, 1:max_age2]
      CAA_hist2[, max_age2] <- rowSums(CAA_hist[, max_age2:max_age, drop = FALSE])
    }

    min_age <- 1 # Increase min-age until no zero's in terminal year
    while(CAA_hist2[nrow(CAA_hist), min_age] <= 0) min_age <- min_age + 1
    if(colSums(CAA_hist2)[1] < 0.01 * sum(CAA_hist2)) min_age <- min_age + 1
    if(min_age > 1) {
      CAA_hist2 <- CAA_hist2[, -c(1:(min_age-1)), drop = FALSE]
      CAA_hist2 <- rowSums(CAA_hist[, 1:min_age, drop = FALSE])
    }
    if(ncol(CAA_hist2) == 1) stop("Only one age class left after consolidating plus- and minus- groups to remove zeros.")

    # Any missing zeros
    ages <- min_age:max_age2

  }

  CAA_hist2[is.na(CAA_hist2) | CAA_hist2 < 1e-8] <- 1e-8
  max_ind <- CAA_hist2[, ncol(CAA_hist2) - 1] == 1e-8
  CAA_hist2[max_ind, ncol(CAA_hist2) - 1] <- CAA_hist2[max_ind, ncol(CAA_hist2)]

  Wa2 <- Wa[min_age:max_age2]
  Wa2[length(Wa2)] <- mean(Wa[max_age2:max_age])
  La2 <- La[min_age:max_age2]
  La2[length(Wa2)] <- mean(La[max_age2:max_age])
  mat_age2 <- mat_age[min_age:max_age2]
  mat_age2[length(mat_age2)] <- mean(mat_age[max_age2:max_age])
  mat_age2 <- mat_age2/max(mat_age2)

  LH <- list(LAA = La2, WAA = Wa2, Linf = Linf, K = K, t0 = t0, a = a, b = b, A50 = A50, A95 = A95)

  data <- list(model = "VPA", I_hist = I_hist, CAA_hist = CAA_hist2, n_y = length(Year), max_age = length(ages),
               M = M[min_age:max_age2], weight = Wa2, mat = mat_age2, vul_type_term = vulnerability,
               I_type = I_type, nitF = as.integer(nitF),
               n_vulpen = vul_pen[1], sigma_vulpen = vul_pen[2], n_Rpen = R_pen[1], sigma_Rpen = R_pen[2])

  # Starting values
  params <- list()
  if(!is.null(start)) {
    if(!is.null(start$F_term) && is.numeric(start$F_term)) params$logF_term <- log(start$F_term)
    if(!is.null(start$F_ratio) && is.numeric(start$F_ratio)) params$logF_ratio <- log(start$F_ratio)
    if(!is.null(start$vul_par) && is.numeric(start$vul_par)) {
      if(vulnerability == "logistic") {
        if(start$vul_par[1] > 0.75 * max_age) stop("start$vul_par[1] needs to be less than 0.75 * Data@MaxAge (see help).")
        if(length(start$vul_par) < 2) stop("Two parameters needed for start$vul_par with logistic vulnerability (see help).")
        if(start$vul_par[1] <= start$vul_par[2]) stop("start$vul_par[1] needs to be greater than start$vul_par[2] (see help).")

        params$vul_par_sep <- c(logit(start$vul_par[1]/max_age/0.75), log(start$vul_par[1] - start$vul_par[2]))
      }
      if(vulnerability == "dome") {
        if(start$vul_par[1] > 0.75 * max_age) stop("start$vul_par[1] needs to be less than 0.75 * Data@MaxAge (see help).")
        if(length(start$vul_par) < 4) stop("Four parameters needed for start$vul_par with dome vulnerability (see help).")
        if(start$vul_par[1] <= start$vul_par[2]) stop("start$vul_par[1] needs to be greater than start$vul_par[2] (see help).")
        if(start$vul_par[3] <= start$vul_par[1] || start$vul_par[3] >= max_age) {
          stop("start$vul_par[3] needs to be between start$vul_par[1] and Data@MaxAge (see help).")
        }
        if(start$vul_par[4] <= 0 || start$vul_par[4] >= 1) stop("start$vul_par[4] needs to be between 0-1 (see help).")

        params$vul_par_sep <- c(logit(start$vul_par[1]/max_age/0.75), log(start$vul_par[1] - start$vul_par[2]),
                                logit(1/(max_age - start$vul_par[1])), logit(start$vul_par[4]))
      }
      if(vulnerability == "free") {
        if(length(start$vul_par) < length(ages)) stop(paste0("start$vul_par needs to be of length", length(ages), "."))
        params$vul_par_sep <- log(start$vul_par[1:length(ages)])
      }
    }
    if(!is.null(start$sigma) && is.numeric(start$sigma)) params$log_sigma <- log(start$sigma)
  }

  if(is.null(params$logF_term)) params$logF_term <- log(0.2)
  if(is.null(params$logF_ratio)) params$logF_ratio <- log(1)
  if(is.null(params$vul_par)) {
    if(vulnerability == "free") {
      params$vul_par <- rep(0, length(ages) - 1)
    } else {
      CAA_mode <- which.max(colSums(CAA_hist, na.rm = TRUE))
      CAA_mode <- ifelse(CAA_mode + 1 > max_age2, max_age2 - 1, CAA_mode)
      if((is.na(Data@LFC[x]) && is.na(Data@LFS[x])) || (Data@LFC[x] > Linf) || (Data@LFS[x] > Linf)) {
        if(vulnerability == "logistic") params$vul_par <- c(logit(CAA_mode/(max_age2 - min_age)/0.75), log(1))
        if(vulnerability == "dome") {
          params$vul_par <- c(logit(CAA_mode/(max_age2 - min_age)/0.75), log(1), logit(1/(max_age2 - CAA_mode)), logit(0.5))
        }
      } else {
        A5 <- min(iVB(t0, K, Linf, Data@LFC[x]), CAA_mode-1)
        Afull <- min(iVB(t0, K, Linf, Data@LFS[x]), 0.5 * max_age2)
        A5 <- min(A5, Afull - 0.5)
        A50_vul <- mean(c(A5, Afull))

        if(vulnerability == "logistic") params$vul_par <- c(logit(Afull/(max_age2 - min_age)/0.75), log(Afull - A50_vul))
        if(vulnerability == "dome") {
          params$vul_par <- c(logit(Afull/max_age2/0.75), log(Afull - A50_vul), logit(1/(max_age2 - Afull)), logit(0.5))
        }
      }
    }
  }
  if(is.null(params$log_sigma)) {
    sigmaI <- max(0.05, sdconv(1, Data@CV_Ind[x]), na.rm = TRUE)
    params$log_sigma <- log(sigmaI)
  }

  info <- list(Year = Year, data = data, params = params, LH = LH, SR = SR, ages = ages, control = control, fix_h = fix_h)

  map <- list()
  if(fix_sigma) map$log_sigma <- factor(NA)
  if(vulnerability == "free") {
    fixed_ind <- round(0.5 * length(ages))
    free_vul_par <- rep(log(1), length(ages) - 1)
    free_vul_par[fixed_ind] <- NA
    free_vul_par[!is.na(free_vul_par)] <- 1:(length(ages)-2)
    map$vul_par <- factor(free_vul_par)
  }
  if(fix_Fratio) map$logF_ratio <- factor(NA)
  if(any(names(dots) == "fix_F")) {
    if(dots$fixF) map$logF_term <- factor(NA)
  }

  obj <- MakeADFun(data = info$data, parameters = info$params, hessian = TRUE,
                   map = map, DLL = "MSEtool", silent = silent)

  mod <- optimize_TMB_model(obj, control, opt_hess, n_restart)
  opt <- mod[[1]]
  SD <- mod[[2]]
  report <- obj$report(obj$env$last.par.best)

  # Terminal-year + 1 abundance
  report <- projection_VPA_internal(report, info, R_pen[1])
  Yearplusone <- c(Year, max(Year) + 1)

  Assessment <- new("Assessment", Model = "VPA",
                    Name = Data@Name, conv = !is.character(SD) && SD$pdHess,
                    FMort = structure(apply(report$F, 1, max), names = Year),
                    B = structure(report$B, names = Yearplusone), SSB = structure(report$E, names = Yearplusone),
                    VB = structure(report$VB, names = Yearplusone),
                    R = structure(report$N[, 1], names = Yearplusone), N = structure(rowSums(report$N), names = Yearplusone),
                    N_at_age = report$N, Selectivity = report$vul, h = NaN,
                    Obs_Catch = structure(if(expanded) colSums(t(CAA_hist2) * data$weight) else C_hist, names = Year),
                    Obs_Index = structure(I_hist, names = Year), Obs_C_at_age = CAA_hist2,
                    Catch = structure(colSums(t(report$CAApred) * data$weight), names = Year),
                    Index = structure(I_hist, names = Year), C_at_age = report$CAApred,
                    NLL = structure(c(report$fn, report$nll_comp[1], sum(report$nll_comp[2:3], report$penalty, report$prior)),
                                    names = c("Total", "Index", "Penalty")),
                    info = info, obj = obj, opt = opt, SD = SD, TMB_report = report,
                    dependencies = dependencies)

  if(Assessment@conv) {

    report$vul_refpt <- apply(report$vul[(length(Year)-vul_pen[1]+1):length(Year), , drop = FALSE], 2, mean)
    report$vul_refpt <- report$vul_refpt/max(report$vul_refpt)

    E <- report$E[1:(length(report$E)-min_age)]
    R <- report$N[(min_age + 1):length(report$E), 1]
    ref_pt <- SCA_refpt_calc(E, R, data$weight, data$mat, data$M, report$vul_refpt, SR, fix_h, ifelse(fix_h, Data@steep[x], NA))

    report <- c(report, ref_pt)

    Z_mat <- t(report$F) + data$M
    VB_mid <- t(report$N[-ncol(report$N), ]) * data$weight * (1 - exp(-Z_mat))/Z_mat

    Assessment@FMSY <- report$FMSY
    Assessment@U <- structure(Assessment@Catch/colSums(VB_mid), names = Year)
    Assessment@UMSY <- report$MSY/report$VBMSY
    Assessment@U_UMSY <- Assessment@U/Assessment@UMSY
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
    Assessment@F_FMSY <- Assessment@FMort/report$FMSY
    Assessment@B_BMSY <- structure(report$B/report$BMSY, names = Yearplusone)
    Assessment@B_B0 <- structure(report$B/report$B0, names = Yearplusone)
    Assessment@SSB_SSBMSY <- structure(report$E/report$EMSY, names = Yearplusone)
    Assessment@SSB_SSB0 <- structure(report$E/report$E0, names = Yearplusone)
    Assessment@VB_VBMSY <- structure(report$VB/report$VBMSY, names = Yearplusone)
    Assessment@VB_VB0 <- structure(report$VB/report$VB0, names = Yearplusone)
    Assessment@TMB_report <- report
    Assessment@info <- info
  }

  return(Assessment)
}




projection_VPA_internal <- function(report, info, nR) {
  max_age <- info$data$max_age
  termY <- nrow(report$N)
  N <- numeric(max_age)
  N[2:max_age] <- report$N[termY, 1:(max_age-1)] * exp(-info$data$M[1:(max_age-1)] - report$F[termY, 1:(max_age-1)])
  N[max_age] <- N[max_age] + report$N[termY, max_age] * exp(-info$data$M[max_age] - report$F[termY, max_age])
  N[1] <- exp(sum(log(report$N[(termY - nR + 1):termY, 1]))/nR)

  report$N <- rbind(report$N, N)
  report$E <- c(report$E, sum(N * info$data$mat * info$data$weight))
  report$VB <- c(report$VB, sum(N * report$vul[termY, ] * info$data$weight))
  report$B <- c(report$B, sum(N * info$data$weight))

  return(report)
}
