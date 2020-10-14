
#' Simple Stock Synthesis
#'
#' A simple age-structured model (\link{SCA_Pope}) fitted to a time series of catch going back to unfished conditions.
#' Terminal depletion (ratio of current biomass to unfished biomass) is by default fixed to 0.4.
#' Selectivity is fixed to the maturity ogive,
#' although it can be overridden with the start argument. The sole parameter estimated is R0 (unfished recruitment).
#'
#' @param x A position in the Data object (by default, equal to one for assessments).
#' @param Data An object of class Data
#' @param dep Depletion value to use in the model. Can be an expression that will be evaluated inside the function.
#' @param SR Stock-recruit function (either \code{"BH"} for Beverton-Holt or \code{"Ricker"}).
#' @param rescale A multiplicative factor that rescales the catch in the assessment model, which
#' can improve convergence. By default, \code{"mean1"} scales the catch so that time series mean is 1, otherwise a numeric.
#' Output is re-converted back to original units.
#' @param start Optional named list of starting values. Entries can be expressions that are evaluated in the function:
#' \itemize{
#' \item R0 - unfished recruitment
#' \item vul_par - a length-two vector for the age of 95\% and 50\% fleet selectivity. Fixed to maturity otherwise.
#' }
#' @param silent Logical, passed to \code{\link[TMB]{MakeADFun}}, whether TMB
#' will print trace information during optimization. Used for dignostics for model convergence.
#' @param opt_hess Logical, whether the hessian function will be passed to \code{\link[stats]{nlminb}} during optimization
#' (this generally reduces the number of iterations to convergence, but is memory and time intensive and does not guarantee an increase
#' in convergence rate).
#' @param n_restart The number of restarts (calls to \code{\link[stats]{nlminb}}) in the optimization procedure, so long as the model
#' hasn't converged. The optimization continues from the parameters from the previous (re)start.
#' @param control A named list of agruments for optimization to be passed to \code{\link[stats]{nlminb}}.
#' @param ... Other arguments to be passed (not currently used).
#' @references
#' Cope, J.M. 2013. Implementing a statistical catch-at-age model (Stock Synthesis) as a tool for
#' deriving overfishing limits in data-limited situations. Fisheries Research 142:3-14.
#' @author Q. Huynh
#' @examples
#' res <- SSS(1, Data = Red_snapper)
#'
#' SSS_MP <- make_MP(SSS, HCR40_10, dep = 0.3) # Always assume depletion = 0.3
#' @useDynLib MSEtool
#' @export
SSS <- function(x = 1, Data, dep = 0.4, SR = c("BH", "Ricker"), rescale = "mean1",
                start = NULL, silent = TRUE, opt_hess = FALSE, n_restart = ifelse(opt_hess, 0, 1),
                control = list(iter.max = 2e5, eval.max = 4e5), ...) {

  dependencies <- "Data@Cat, Data@steep, Data@Mort, Data@L50, Data@L95, Data@vbK, Data@vbLinf, Data@vbt0, Data@wla, Data@wlb, Data@MaxAge"
  dots <- list(...)
  start <- lapply(start, eval, envir = environment())
  dep <- eval(dep)
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

  n_y <- length(C_hist)
  I_hist <- rep(NA_real_, n_y)
  I_hist[1] <- 1
  I_hist[n_y] <- dep

  max_age <- as.integer(-log(0.01)/Data@Mort[x])
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

  if(rescale == "mean1") rescale <- 1/mean(C_hist)
  data <- list(model = "SCA_Pope", C_hist = C_hist, rescale = rescale, I_hist = I_hist,
               CAA_hist = matrix(0, n_y, max_age), CAA_n = rep(0, n_y), n_y = n_y, max_age = max_age, M = M,
               weight = Wa, mat = mat_age, vul_type = "logistic", I_type = "B",
               SR_type = SR, CAA_dist = "multinomial", est_early_rec_dev = rep(0, max_age - 1),
               est_rec_dev = rep(0, n_y))

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
    if(!is.null(start$vul_par) && is.numeric(start$vul_par)) {
      if(start$vul_par[1] > 0.75 * max_age) stop("start$vul_par[1] needs to be less than 0.75 * Data@MaxAge (see help).")
      if(length(start$vul_par) < 2) stop("Two parameters needed for start$vul_par with logistic vulnerability (see help).")
      if(start$vul_par[1] <= start$vul_par[2]) stop("start$vul_par[1] needs to be greater than start$vul_par[2] (see help).")

      params$vul_par <- c(logit(start$vul_par[1]/max_age/0.75), log(start$vul_par[1] - start$vul_par[2]))
    }
  }

  if(is.null(params$R0x)) {
    params$R0x <- ifelse(is.null(Data@OM$R0[x]), log(mean(data$C_hist)) + 4, log(1.5 * rescale * Data@OM$R0[x]))
  }
  if(is.null(params$transformed_h)) {
    h_start <- Data@steep[x]
    if(SR == "BH") {
      h_start <- (h_start - 0.2)/0.8
      params$transformed_h <- logit(h_start)
    } else {
      params$transformed_h <- log(h_start - 0.2)
    }
  }
  if(is.null(params$vul_par)) params$vul_par <- c(logit(min(A95, 0.74 * max_age)/max_age/0.75), log(A95-A50))

  params$U_equilibrium <- 0
  params$log_sigma <- params$log_tau <- log(0.01)
  params$log_early_rec_dev <- rep(0, max_age - 1)
  params$log_rec_dev <- rep(0, n_y)

  info <- list(Year = Year, data = data, params = params, LH = LH, control = control)

  map <- list()
  map$transformed_h <- map$U_equilibrium <- map$log_sigma <- map$log_tau <- factor(NA)
  map$vul_par <- factor(c(NA, NA))
  map$log_early_rec_dev <- factor(rep(NA, max_age - 1))
  map$log_rec_dev <- factor(rep(NA, n_y))

  obj <- MakeADFun(data = info$data, parameters = info$params, hessian = TRUE,
                   map = map, DLL = "MSEtool", silent = silent)

  # Add starting values for rec-devs and increase R0 start value if U is too high (> 0.975)
  high_U <- try(obj$report(obj$par)$penalty > 0, silent = TRUE)
  if(!is.character(high_U) && high_U) {
    while(obj$par["R0x"] < 30 && obj$report(obj$par)$penalty > 0) {
      obj$par["R0x"] <- obj$par["R0x"] + 1
    }
  }

  mod <- optimize_TMB_model(obj, control, opt_hess, n_restart)
  opt <- mod[[1]]
  SD <- mod[[2]]
  report <- obj$report(obj$env$last.par.best)

  Yearplusone <- c(Year, max(Year) + 1)

  nll_report <- ifelse(is.character(opt), ifelse(integrate, NA, report$nll), opt$objective)
  Assessment <- new("Assessment", Model = "SSS", Name = Data@Name, conv = !is.character(SD) && SD$pdHess,
                    B0 = report$B0, R0 = report$R0, N0 = report$N0,
                    SSB0 = report$E0, VB0 = report$VB0,
                    h = report$h, U = structure(report$U, names = Year),
                    B = structure(report$B, names = Yearplusone),
                    B_B0 = structure(report$B/report$B0, names = Yearplusone),
                    SSB = structure(report$E, names = Yearplusone),
                    SSB_SSB0 = structure(report$E/report$E0, names = Yearplusone),
                    VB = structure(report$VB, names = Yearplusone),
                    VB_VB0 = structure(report$VB/report$VB0, names = Yearplusone),
                    R = structure(report$R, names = Yearplusone),
                    N = structure(rowSums(report$N), names = Yearplusone),
                    N_at_age = report$N,
                    Selectivity = matrix(report$vul, nrow = length(Year),
                                         ncol = max_age, byrow = TRUE),
                    Obs_Catch = structure(C_hist, names = Year),
                    Obs_Index = structure(I_hist, names = Year),
                    Catch = structure(colSums(t(report$CAApred) * Wa), names = Year),
                    Index = structure(report$Ipred, names = Year),
                    C_at_age = report$CAApred,
                    NLL = structure(nll_report, names = "Total"),
                    info = info, obj = obj, opt = opt, SD = SD, TMB_report = report,
                    dependencies = dependencies)

  if(Assessment@conv) {
    ref_pt <- SCA_Pope_MSY_calc(Arec = report$Arec, Brec = report$Brec, M = M, weight = Wa, mat = mat_age, vul = report$vul, SR = SR)
    report <- c(report, ref_pt)

    Assessment@UMSY <- report$UMSY
    Assessment@MSY <- report$MSY
    Assessment@BMSY <- report$BMSY
    Assessment@SSBMSY <- report$EMSY
    Assessment@VBMSY <- report$VBMSY
    Assessment@U_UMSY <- structure(report$U/report$UMSY, names = Year)
    Assessment@B_BMSY <- structure(report$B/report$BMSY, names = Yearplusone)
    Assessment@SSB_SSBMSY <- structure(report$E/report$EMSY, names = Yearplusone)
    Assessment@VB_VBMSY <- structure(report$VB/report$VBMSY, names = Yearplusone)
    Assessment@TMB_report <- report
  }
  return(Assessment)
}
class(SSS) <- "Assess"


