
#' Simple Stock Synthesis
#'
#' A simple age-structured model ([SCA_Pope]) fitted to a time series of catch going back to unfished conditions.
#' Terminal depletion (ratio of current total biomass to unfished biomass) is by default fixed to 0.4. Selectivity is fixed 
#' to the maturity ogive, although it can be overridden with the start argument. The sole parameter estimated is 
#' R0 (unfished recruitment), with no process error.
#'
#' @param x A position in the Data object (by default, equal to one for assessments).
#' @param Data An object of class Data
#' @param dep Depletion value to use in the model. Can be an expression that will be evaluated inside the function.
#' @param SR Stock-recruit function (either `"BH"` for Beverton-Holt or `"Ricker"`).
#' @param rescale A multiplicative factor that rescales the catch in the assessment model, which
#' can improve convergence. By default, `"mean1"` scales the catch so that time series mean is 1, otherwise a numeric.
#' Output is re-converted back to original units.
#' @param start Optional named list of starting values. Entries can be expressions that are evaluated in the function:
#' \itemize{
#' \item `R0` Unfished recruitment
#' \item `vul_par` A length-two vector for the age of 95% and 50% fleet selectivity. Fixed to maturity otherwise.
#' }
#' @param prior A named list for the parameters of any priors to be added to the model. See details in `SCA_Pope`.
#' @param silent Logical, passed to [TMB::MakeADFun()], whether TMB
#' will print trace information during optimization. Used for diagnostics for model convergence.
#' @param opt_hess Logical, whether the hessian function will be passed to [stats::nlminb()] during optimization
#' (this generally reduces the number of iterations to convergence, but is memory and time intensive and does not guarantee an increase
#' in convergence rate).
#' @param n_restart The number of restarts (calls to [stats::nlminb()]) in the optimization procedure, so long as the model
#' hasn't converged. The optimization continues from the parameters from the previous (re)start.
#' @param control A named list of arguments for optimization to be passed to [stats::nlminb()].
#' @param ... Other arguments to be passed (not currently used).
#' @details In SAMtool, SSS is an implementation of [SCA_Pope] with fixed final depletion 
#' (in terms of total biomass, not spawning biomass) assumption.
#' @references
#' Cope, J.M. 2013. Implementing a statistical catch-at-age model (Stock Synthesis) as a tool for
#' deriving overfishing limits in data-limited situations. Fisheries Research 142:3-14.
#' @author Q. Huynh
#' @return An object of class [Assessment-class].
#' @examples
#' res <- SSS(Data = Red_snapper)
#'
#' SSS_MP <- make_MP(SSS, HCR40_10, dep = 0.3) # Always assume depletion = 0.3
#' @useDynLib SAMtool
#' @export
SSS <- function(x = 1, Data, dep = 0.4, SR = c("BH", "Ricker"), 
                rescale = "mean1", start = NULL, prior = list(),
                silent = TRUE, opt_hess = FALSE, n_restart = ifelse(opt_hess, 0, 1),
                control = list(iter.max = 2e5, eval.max = 4e5), ...) {
  
  catch_eq <- "Pope"
  dependencies <- "Data@Cat, Data@steep, Data@Mort, Data@L50, Data@L95, Data@LFC, Data@LFS, Data@vbK, Data@vbLinf, Data@vbt0, Data@wla, Data@wlb, Data@MaxAge"
  dots <- list(...)
  start <- lapply(start, eval, envir = environment())
  dep <- eval(dep)
  SR <- match.arg(SR)

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
  I_hist <- matrix(NA_real_, n_y, 1)
  I_hist[1] <- 1
  I_hist[n_y] <- dep

  max_age <- as.integer(-log(0.01)/Data@Mort[x])
  n_age <- max_age + 1
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
  Wa <- a * La ^ b
  A50 <- min(0.5 * max_age, iVB(t0, K, Linf, Data@L50[x]))
  A95 <- max(A50+0.5, iVB(t0, K, Linf, Data@L95[x]))
  mat_age <- c(0, 1/(1 + exp(-log(19) * (c(1:max_age) - A50)/(A95 - A50))))
  mat_age <- mat_age/max(mat_age)
  LH <- list(LAA = La, WAA = Wa, Linf = Linf, K = K, t0 = t0, a = a, b = b, A50 = A50, A95 = A95)
  
  # Generate priors
  prior <- make_prior(prior, nsurvey = 0, ifelse(SR == "BH", 1, 2), msg = FALSE)
  if (rescale == "mean1") rescale <- 1/mean(C_hist)
  
  data <- list(model = "SCA", C_hist = C_hist, rescale = rescale, 
               I_hist = I_hist, I_sd = matrix(0.01, n_y, 1), I_units = 1, I_vul = matrix(1, n_age, 1), 
               abs_I = 0, nsurvey = 1, LWT = 1,
               CAA_hist = matrix(0, n_y, n_age), CAA_n = rep(0, n_y), 
               CAL_hist = matrix(0, n_y, 1), CAL_n = rep(0, n_y),
               IAA_hist = array(0, c(n_y, n_age, 1)), IAA_n = matrix(0, n_y, 1),
               n_y = n_y, n_age = n_age, n_bin = 1L, M_data = 1,
               weight = Wa, PLA = matrix(1, n_age, 1), mat = mat_age, vul_type = "logistic",
               SR_type = SR, comp_dist = "multinomial", catch_eq = catch_eq,
               est_early_rec_dev = rep(0, n_age - 1), est_rec_dev = rep(0, n_y), yindF = 0,
               tv_M = "none", M_bounds = c(0, 1e4), use_prior = prior$use_prior, prior_dist = prior$pr_matrix,
               sim_process_error = 0L)
  if (any(names(dots) == "M_at_age") && dots$M_at_age) data$M_data <- M

  # Starting values
  params <- list()
  if (!is.null(start)) {
    if (!is.null(start$R0) && is.numeric(start$R0)) params$R0x <- log(start$R0[1] * rescale)
    if (!is.null(start$h) && is.numeric(start$h)) {
      if (SR == "BH") {
        h_start <- (start$h[1] - 0.2)/0.8
        params$transformed_h <- logit(h_start)
      } else {
        params$transformed_h <- log(start$h[1] - 0.2)
      }
    }
    if (!is.null(start$vul_par) && is.numeric(start$vul_par)) {
      if (start$vul_par[1] > 0.75 * max_age) stop("start$vul_par[1] needs to be less than 0.75 * Data@MaxAge (see help).")
      if (length(start$vul_par) < 2) stop("Two parameters needed for start$vul_par with logistic vulnerability (see help).")
      if (start$vul_par[1] <= start$vul_par[2]) stop("start$vul_par[1] needs to be greater than start$vul_par[2] (see help).")

      params$vul_par <- c(logit(start$vul_par[1]/max_age/0.75), log(start$vul_par[1] - start$vul_par[2]))
    }
  }

  if (is.null(params$R0x)) {
    params$R0x <- ifelse(is.null(Data@OM$R0[x]), log(mean(data$C_hist)) + 4, log(1.5 * rescale * Data@OM$R0[x]))
  }
  if (is.null(params$transformed_h)) {
    h_start <- Data@steep[x]
    if (SR == "BH") {
      h_start <- (h_start - 0.2)/0.8
      params$transformed_h <- logit(h_start)
    } else {
      params$transformed_h <- log(h_start - 0.2)
    }
  }
  params$log_M0 <- log(M) %>% mean()
  params$logit_M_walk <- rep(0, n_y)
  params$F_equilibrium <- 0 
  
  if (is.null(params$vul_par)) params$vul_par <- c(logit(min(A95, 0.74 * max_age)/max_age/0.75), log(A95-A50))

  params$log_F_dev <- rep(0, n_y)
  params$log_omega <- params$log_tau <- params$log_tau_M <- log(0.01)
  params$log_early_rec_dev <- rep(0, n_age - 1)
  params$log_rec_dev <- rep(0, n_y)

  info <- list(Year = Year, data = data, params = params, LH = LH, control = control)

  map <- list()
  if (!prior$use_prior[2]) map$transformed_h <- factor(NA)
  if (!prior$use_prior[3]) map$log_M0 <- factor(NA)
  map$logit_M_walk <- factor(rep(NA, n_y))
  map$F_equilibrium <- factor(NA)
  map$vul_par <- factor(c(NA, NA))
  map$log_F_dev <- factor(rep(NA, n_y))
  map$log_omega <- map$log_tau <- map$log_tau_M <- factor(NA)
  map$log_early_rec_dev <- factor(rep(NA, n_age - 1))
  map$log_rec_dev <- factor(rep(NA, n_y))

  obj <- MakeADFun(data = info$data, parameters = info$params, hessian = TRUE,
                   map = map, DLL = "SAMtool", silent = silent)

  # Add starting values for rec-devs and increase R0 start value if U is too high (> 0.975)
  high_U <- try(obj$report(obj$par)$penalty > 0, silent = TRUE)
  if (!is.character(high_U) && !is.na(high_U) && high_U) {
    while (obj$par["R0x"] < 30 && obj$report(obj$par)$penalty > 0) {
      obj$par["R0x"] <- obj$par["R0x"] + 1
    }
  }

  mod <- optimize_TMB_model(obj, control, opt_hess, n_restart)
  opt <- mod[[1]]
  SD <- mod[[2]]
  report <- obj$report(obj$env$last.par.best)

  Yearplusone <- c(Year, max(Year) + 1)
  
  report$dynamic_SSB0 <- SCA_dynamic_SSB0(obj) %>% 
    structure(names = Yearplusone)
  
  nll_report <- ifelse(is.character(opt), ifelse(integrate, NA, report$nll), opt$objective)
  Assessment <- new("Assessment", Model = "SSS", 
                    Name = Data@Name, conv = SD$pdHess,
                    B0 = report$B0, R0 = report$R0, N0 = report$N0,
                    SSB0 = report$E0, VB0 = report$VB0,
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
                                         ncol = n_age, byrow = TRUE),
                    Obs_Catch = structure(C_hist, names = Year),
                    Obs_Index = structure(I_hist, dimnames = list(Year, "Depletion")),
                    Catch = structure(report$Cpred, names = Year),
                    Index = structure(report$Ipred, dimnames = list(Year, "Depletion")),
                    C_at_age = report$CAApred,
                    NLL = structure(nll_report, names = "Total"),
                    info = info, obj = obj, opt = opt, SD = SD, TMB_report = report,
                    dependencies = dependencies)
  Assessment@h <- report$h
  Assessment@U <- structure(report$U, names = Year)
  
  if (Assessment@conv) {
    ref_pt <- ref_pt_SCA(obj = obj, report = report)
    report <- c(report, ref_pt[1:6])
    
    Assessment@UMSY <- report$UMSY
    Assessment@U_UMSY <- structure(report$U/report$UMSY, names = Year)
    Assessment@MSY <- report$MSY
    Assessment@BMSY <- report$BMSY
    Assessment@SSBMSY <- report$EMSY
    Assessment@VBMSY <- report$VBMSY
    Assessment@B_BMSY <- structure(report$B/report$BMSY, names = Yearplusone)
    Assessment@SSB_SSBMSY <- structure(report$E/report$EMSY, names = Yearplusone)
    Assessment@VB_VBMSY <- structure(report$VB/report$VBMSY, names = Yearplusone)
    Assessment@TMB_report <- report
    
    catch_eq_fn <- function(Ftarget) {
      projection_SCA(Assessment, Ftarget = Ftarget, p_years = 1, p_sim = 1, obs_error = list(matrix(1, 1, 1), matrix(1, 1, 1)),
                     process_error = matrix(1, 1, 1)) %>% slot("Catch") %>% as.vector()
    }
    Assessment@forecast <- list(per_recruit = ref_pt[[7]], catch_eq = catch_eq_fn)
  }
  return(Assessment)
}
class(SSS) <- "Assess"


