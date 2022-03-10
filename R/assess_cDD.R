#' Continuous Delay-differential assessment model
#'
#' A catch and index-based assessment model. Compared to the discrete delay-difference (annual time-step in production and fishing), the
#' delay-differential model (cDD) is based on continuous recruitment and fishing mortality within a time-step. The continuous model works
#' much better for populations with high turnover (e.g. high F or M, continuous reproduction). This model is conditioned on catch and fits
#' to the observed index. In the state-space version (cDD_SS), recruitment deviations from the stock-recruit relationship are estimated.
#'
#' @param x An index for the objects in \code{Data} when running in closed loop simulation.
#' Otherwise, equals to 1 when running an assessment.
#' @param Data An object of class \linkS4class{Data}.
#' @param AddInd A vector of integers or character strings indicating the indices to be used in the model. Integers assign the index to
#' the corresponding index in Data@@AddInd, "B" (or 0) represents total biomass in Data@@Ind, "VB" represents vulnerable biomass in
#' Data@@VInd, and "SSB" represents spawning stock biomass in Data@@SpInd.
#' @param SR Stock-recruit function (either \code{"BH"} for Beverton-Holt or \code{"Ricker"}).
#' @param rescale A multiplicative factor that rescales the catch in the assessment model, which
#' can improve convergence. By default, \code{"mean1"} scales the catch so that time series mean is 1, otherwise a numeric.
#' Output is re-converted back to original units.
#' @param MW Logical, whether to fit to mean weight. In closed-loop simulation, mean weight will be grabbed from \code{Data@@Misc[[x]]$MW},
#' otherwise calculated from \code{Data@@CAL}.
#' @param start Optional list of starting values. Entries can be expressions that are evaluated in the function. See details.
#' @param prior A named list for the parameters of any priors to be added to the model. See below.
#' @param fix_h Logical, whether to fix steepness to value in \code{Data@@steep} in the assessment model.
#' @param fix_sigma Logical, whether the standard deviation of the index is fixed. If \code{TRUE},
#' sigma is fixed to value provided in \code{start} (if provided), otherwise, value based on \code{Data@@CV_Ind}.
#' @param fix_tau Logical, the standard deviation of the recruitment deviations is fixed. If \code{TRUE},
#' tau is fixed to value provided in \code{start} (if provided), otherwise, equal to 1.
#' @param dep The initial depletion in the first year of the model. A tight prior is placed on the model objective function
#' to estimate the equilibrium fishing mortality corresponding to the initial depletion. Due to this tight prior, this F
#' should not be considered to be an independent model parameter. Set to zero to eliminate this prior.
#' @param integrate Logical, whether the likelihood of the model integrates over the likelihood
#' of the recruitment deviations (thus, treating it as a state-space variable). Otherwise, recruitment deviations are penalized parameters.
#' @param LWT A named list of likelihood weights. For \code{LWT$Index}, a vector of likelihood weights for each survey, while
#' for \code{LWT$MW} a numeric.
#' @param silent Logical, passed to \code{\link[TMB]{MakeADFun}}, whether TMB
#' will print trace information during optimization. Used for diagnostics for model convergence.
#' @param n_itF Integer, the number of iterations to solve F conditional on the observed catch.
#' @param opt_hess Logical, whether the hessian function will be passed to \code{\link[stats]{nlminb}} during optimization
#' (this generally reduces the number of iterations to convergence, but is memory and time intensive and does not guarantee an increase
#' in convergence rate). Ignored if \code{integrate = TRUE}.
#' @param n_restart The number of restarts (calls to \code{\link[stats]{nlminb}}) in the optimization procedure, so long as the model
#' hasn't converged. The optimization continues from the parameters from the previous (re)start.
#' @param control A named list of parameters regarding optimization to be passed to
#' \code{\link[stats]{nlminb}}.
#' @param inner.control A named list of arguments for optimization of the random effects, which
#' is passed on to \code{\link[TMB]{newton}} via \code{\link[TMB]{MakeADFun}}.
#' @param ... Additional arguments (not currently used).
#' @return An object of \code{\linkS4class{Assessment}} containing objects and output
#' from TMB.
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
#' For \code{start} (optional), a named list of starting values of estimates can be provided for:
#' \itemize{
#' \item \code{R0} Unfished recruitment. Otherwise, Data@@OM$R0[x] is used in closed-loop, and 400\% of mean catch otherwise.
#' \item \code{h} Steepness. Otherwise, Data@@steep[x] is used, or 0.9 if empty.
#' \item \code{Kappa} Delay-differential Kappa parameter. Otherwise, calculated from biological parameters in the Data object.
#' \item \code{F_equilibrium} Equilibrium fishing mortality leading into first year of the model (to determine initial depletion). By default, 0.
#' \item \code{tau} Lognormal SD of the recruitment deviations (process error) for \code{DD_SS}. By default, Data@@sigmaR[x].
#' \item \code{sigma} Lognormal SD of the index (observation error). By default, Data@@CV_Ind[x]. Not
#' used if multiple indices are used.
#' \item \code{sigma_W} Lognormal SD of the mean weight (observation error). By default, 0.1.
#' }
#' 
#' Multiple indices are supported in the model. Data@@Ind, Data@@VInd, and Data@@SpInd are all assumed to be biomass-based.
#' For Data@@AddInd, Data@@I_units are used to identify a biomass vs. abundance-based index.
#'
#' @section Online Documentation:
#' Model description and equations are available on the openMSE 
#' \href{https://openmse.com/features-assessment-models/1-dd/}{website}.
#' 
#' @author Q. Huynh
#' @references
#' Hilborn, R., and Walters, C., 1992. Quantitative Fisheries Stock Assessment: Choice,
#' Dynamics and Uncertainty. Chapman and Hall, New York.
#' @section Required Data:
#' \itemize{
#' \item \code{cDD}: Cat, Ind, Mort, L50, vbK, vbLinf, vbt0, wla, wlb, MaxAge
#' \item \code{cDD_SS}: Cat, Ind, Mort, L50, vbK, vbLinf, vbt0, wla, wlb, MaxAge
#' }
#' @section Optional Data:
#' \itemize{
#' \item \code{cDD}: steep
#' \item \code{cDD_SS}: steep, CV_Ind, sigmaR
#' }
#' @examples
#' #### Observation-error delay difference model
#' res <- cDD(Data = MSEtool::Red_snapper)
#'
#' ### State-space version
#' ### Also set recruitment variability SD = 0.6 (since fix_tau = TRUE)
#' res <- cDD_SS(Data = MSEtool::Red_snapper, start = list(tau = 0.6))
#'
#' summary(res@@SD) # Parameter estimates
#' @seealso \link{DD_TMB} \link{plot.Assessment} \link{summary.Assessment} \link{retrospective} \link{profile} \link{make_MP}
#' @export
cDD <- function(x = 1, Data, AddInd = "B", SR = c("BH", "Ricker"), rescale = "mean1", MW = FALSE, start = NULL, prior = list(), fix_h = TRUE,
                dep = 1, LWT = list(), n_itF = 5L, silent = TRUE, opt_hess = FALSE, n_restart = ifelse(opt_hess, 0, 1),
                control = list(iter.max = 5e3, eval.max = 1e4), ...) {
  cDD_(x = x, Data = Data, AddInd = AddInd, state_space = FALSE, SR = SR, rescale = rescale, MW = MW, start = start, prior = prior,
       fix_h = fix_h, fix_sigma = FALSE, fix_tau = TRUE, dep = dep, LWT = LWT, n_itF = n_itF,
       integrate = FALSE, silent = silent, opt_hess = opt_hess, n_restart = n_restart,
       control = control, inner.control = list(), ...)
}
class(cDD) <- "Assess"


#' @rdname cDD
#' @export
cDD_SS <- function(x = 1, Data, AddInd = "B", SR = c("BH", "Ricker"), rescale = "mean1", MW = FALSE, start = NULL, prior = list(),
                   fix_h = TRUE, fix_sigma = FALSE, fix_tau = TRUE, dep = 1, LWT = list(), n_itF = 5L,
                   integrate = FALSE, silent = TRUE, opt_hess = FALSE, n_restart = ifelse(opt_hess, 0, 1),
                   control = list(iter.max = 5e3, eval.max = 1e4), inner.control = list(), ...) {
  cDD_(x = x, Data = Data, AddInd = AddInd, state_space = TRUE, SR = SR, rescale = rescale, MW = MW, start = start, prior = prior,
       fix_h = fix_h, fix_sigma = fix_sigma, fix_tau = fix_tau, dep = dep, LWT = LWT, n_itF = n_itF,
       integrate = integrate, silent = silent, opt_hess = opt_hess, n_restart = n_restart,
       control = control, inner.control = inner.control, ...)
}
class(cDD_SS) <- "Assess"

#' @useDynLib SAMtool
cDD_ <- function(x = 1, Data, AddInd = "B", state_space = FALSE, SR = c("BH", "Ricker"), rescale = "mean1", MW = FALSE, start = NULL,
                 prior = list(), fix_h = TRUE, fix_sigma = FALSE, fix_tau = TRUE, dep = 1, LWT = list(), n_itF = 5L,
                 integrate = FALSE, silent = TRUE, opt_hess = FALSE, n_restart = ifelse(opt_hess, 0, 1),
                 control = list(iter.max = 5e3, eval.max = 1e4), inner.control = list(), ...) {
  dependencies <- "Data@Cat, Data@Ind, Data@Mort, Data@L50, Data@vbK, Data@vbLinf, Data@vbt0, Data@wla, Data@wlb, Data@MaxAge"
  dots <- list(...)
  start <- lapply(start, eval, envir = environment())

  SR <- match.arg(SR)
  Winf <- Data@wla[x] * Data@vbLinf[x]^Data@wlb[x]
  age <- 1:Data@MaxAge
  la <- Data@vbLinf[x] * (1 - exp(-Data@vbK[x] * ((age - Data@vbt0[x]))))
  wa <- Data@wla[x] * la^Data@wlb[x]
  a50V <- iVB(Data@vbt0[x], Data@vbK[x], Data@vbLinf[x],  Data@L50[x])
  a50V <- max(a50V, 1)
  if(any(names(dots) == "yind")) {
    yind <- eval(dots$yind)
  } else {
    ystart <- which(!is.na(Data@Cat[x, ]))[1]
    yind <- ystart:length(Data@Cat[x, ])
  }
  Year <- Data@Year[yind]
  C_hist <- Data@Cat[x, yind]
  ny <- length(C_hist)
  if(any(is.na(C_hist))) stop('Model is conditioned on complete catch time series, but there is missing catch.')

  Ind <- lapply(AddInd, Assess_I_hist, Data = Data, x = x, yind = yind)
  I_hist <- vapply(Ind, getElement, numeric(ny), "I_hist")
  I_sd <- vapply(Ind, getElement, numeric(ny), "I_sd")
  I_units <- vapply(Ind, getElement, numeric(1), "I_units")
  if(is.null(I_hist)) stop("No indices found.", call. = FALSE)
  nsurvey <- ncol(I_hist)
  
  if(MW) {
    if(!is.null(Data@Misc[[x]]$MW)) {
      MW_hist <- Data@Misc[[x]]$MW
    } else {
      MW_hist <- apply(Data@CAL[x, , ], 1, function(xx) {
        weighted.mean(x = Data@wla[x]*Data@CAL_mids^Data@wlb[x], w = xx, na.rm = TRUE)
      })
    }
    MW_hist[MW_hist <= 0] <- NA_real_
  } else {
    MW_hist <- rep(NA_real_, ny)
  }
  
  # Generate priors
  prior <- make_prior(prior, nsurvey, ifelse(SR == "BH", 1, 2), msg = FALSE)

  k <- ceiling(a50V)  # get age nearest to 50% vulnerability (ascending limb)
  k[k > Data@MaxAge/2] <- ceiling(Data@MaxAge/2)  # to stop stupidly high estimates of age at 50% vulnerability
  wk <- wa[k]

  wt_df <- data.frame(t = age[-c(1:(k-1))], W = wa[-c(1:(k-1))], Winf = Winf)
  wt_df$W2 <- c(wt_df$W[2:nrow(wt_df)], NA)
  wt_df <- wt_df[-nrow(wt_df), ]
  
  mod_formula <- formula(W2 ~ Winf + (W - Winf) * exp(-Kappa))
  fit_mod <- nls(mod_formula, wt_df, start = list(Kappa = Data@vbK[x]))
  
  if(!is.null(start$Kappa)) {
    Kappa <- start$Kappa
  } else {
    Kappa <- coef(fit_mod)[["Kappa"]]
  }
  M <- Data@Mort[x]

  if(!state_space && (nsurvey == 1 & AddInd == "B")) fix_sigma <- FALSE # Override: estimate sigma if there's a single survey
  if(rescale == "mean1") rescale <- 1/mean(C_hist)
  if(dep <= 0 || dep > 1) stop("Initial depletion (dep) must be > 0 and <= 1.")
  if(!is.list(LWT)) {
    if(!is.null(LWT) && length(LWT) != nsurvey) stop("LWT needs to be a vector of length ", nsurvey)
    LWT <- list(Index = LWT)
    LWT$MW <- 1
  } else {
    if(is.null(LWT$Index)) LWT$Index <- rep(1, nsurvey)
    if(is.null(LWT$MW)) LWT$MW <- 1 
  }

  data <- list(model = "cDD", Winf = Winf, Kappa = Kappa, ny = ny, k = k, wk = wk, C_hist = C_hist, dep = dep,
               rescale = rescale, I_hist = I_hist, I_units = I_units, I_sd = I_sd, MW_hist = MW_hist,
               SR_type = SR, n_itF = n_itF, LWT = c(LWT$Index, LWT$MW), nsurvey = nsurvey,
               fix_sigma = as.integer(fix_sigma), state_space = as.integer(state_space),
               use_prior = prior$use_prior, prior_dist = prior$pr_matrix)
  LH <- list(LAA = la, WAA = wa, maxage = Data@MaxAge, A50 = k, fit_mod = fit_mod)

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
    if(!is.null(start$M) && is.numeric(start$M)) params$log_M <- log(start$M[1])
    if(!is.null(start$F_equilibrium) && is.numeric(start$F_equilibrium)) params$F_equilibrium <- start$F_equililbrium
    if(!is.null(start[["sigma"]]) && is.numeric(start[["sigma"]])) params$log_sigma <- log(start[["sigma"]])
    if(!is.null(start[["sigma_W"]]) && is.numeric(start[["sigma_W"]])) params$log_sigma_W <- log(start[["sigma_W"]])
    if(!is.null(start$tau) && is.numeric(start$tau)) params$log_tau <- log(start$tau[1])
  }
  
  if(is.null(params$R0x)) {
    params$R0x <- ifelse(is.null(Data@OM$R0[x]), log(4 * mean(data$C_hist)), log(1.5 * rescale * Data@OM$R0[x]))
  }
  if(is.null(params$transformed_h)) {
    h_start <- ifelse(is.na(Data@steep[x]), 0.9, Data@steep[x])
    if(SR == "BH") {
      h_start <- (h_start - 0.2)/0.8
      params$transformed_h <- logit(h_start)
    } else {
      params$transformed_h <- log(h_start - 0.2)
    }
  }
  if(is.null(params$log_M)) params$log_M <- log(M)
  if(is.null(params$F_equilibrium)) params$F_equilibrium <- ifelse(dep < 1, 0.1, 0)
  if(is.null(params[["log_sigma"]])) params$log_sigma <- max(0.05, sdconv(1, Data@CV_Ind[x]), na.rm = TRUE) %>% log()
  if(is.null(params[["log_sigma_W"]])) params$log_sigma_W <- log(0.1)
  if(is.null(params$log_tau)) params$log_tau <- ifelse(is.na(Data@sigmaR[x]), 0.6, Data@sigmaR[x]) %>% log()
  params$log_rec_dev <- rep(0, ny)

  info <- list(Year = Year, data = data, params = params, LH = LH, control = control, inner.control = inner.control)

  map <- list()
  if(fix_h && !prior$use_prior[2]) map$transformed_h <- factor(NA)
  if(!prior$use_prior[3]) map$log_M <- factor(NA)
  if(dep == 1) map$F_equilibrium <- factor(NA)
  if(fix_sigma) map$log_sigma <- factor(NA)
  map$log_sigma_W <- factor(NA)
  if(fix_tau) map$log_tau <- factor(NA)
  if(!state_space) map$log_rec_dev <- factor(rep(NA, ny))

  random <- NULL
  if(integrate) random <- "log_rec_dev"

  obj <- MakeADFun(data = info$data, parameters = info$params, random = random, map = map, hessian = TRUE,
                   DLL = "SAMtool", inner.control = inner.control, silent = silent)
  
  high_F <- try(obj$report(c(obj$par, obj$env$last.par[obj$env$random]))$penalty > 0 ||
                  any(is.na(obj$report(c(obj$par, obj$env$last.par[obj$env$random]))$F)), silent = TRUE)
  if(!is.character(high_F) && !is.na(high_F) && high_F) {
    for(ii in 1:10) {
      obj$par["R0x"] <- 0.5 + obj$par["R0x"]
      if(all(!is.na(obj$report(obj$par)$F)) && 
         obj$report(c(obj$par, obj$env$last.par[obj$env$random]))$penalty == 0) break
    }
  }

  mod <- optimize_TMB_model(obj, control, opt_hess, n_restart)
  opt <- mod[[1]]
  SD <- mod[[2]]
  report <- obj$report(obj$env$last.par.best)

  Yearplusone <- c(Year, max(Year) + 1)
  Yearplusk <- c(Year, max(Year) + 1:k)

  nll_report <- ifelse(is.character(opt), ifelse(integrate, NA, report$nll), opt$objective)
  
  report$dynamic_SSB0 <- cDD_dynamic_SSB0(obj, data = info$data, params = info$params, map = map) %>% 
    structure(names = Yearplusone)
  
  Assessment <- new("Assessment", Model = ifelse(state_space, "cDD_SS", "cDD"),
                    Name = Data@Name, conv = SD$pdHess,
                    B0 = report$B0, R0 = report$R0, N0 = report$N0,
                    SSB0 = report$B0, VB0 = report$B0, h = report$h,
                    FMort = structure(report$F, names = Year),
                    B = structure(report$B, names = Yearplusone),
                    B_B0 = structure(report$B/report$B0, names = Yearplusone),
                    SSB = structure(report$B, names = Yearplusone),
                    SSB_SSB0 = structure(report$B/report$B0, names = Yearplusone),
                    VB = structure(report$B, names = Yearplusone),
                    VB_VB0 = structure(report$B/report$B0, names = Yearplusone),
                    R = structure(report$R, names = Yearplusk),
                    N = structure(report$N, names = Yearplusone),
                    Obs_Catch = structure(C_hist, names = Year),
                    Obs_Index = structure(I_hist, dimnames = list(Year, paste0("Index_", 1:nsurvey))),
                    Catch = structure(report$Cpred, names = Year),
                    Index = structure(report$Ipred, dimnames = list(Year, paste0("Index_", 1:nsurvey))),
                    NLL = structure(c(nll_report, report$nll_comp, report$prior, report$penalty),
                                    names = c("Total", paste0("Index_", 1:nsurvey), "Dev", "Prior", "Penalty")),
                    info = info, obj = obj, opt = opt, SD = SD, TMB_report = report,
                    dependencies = dependencies)

  if(state_space) {
    Assessment@Dev <- structure(report$log_rec_dev, names = Year)
    Assessment@Dev_type <- "log-Recruitment deviations"
  }

  if(Assessment@conv) {
    ref_pt <- ref_pt_cDD(info$data, report$Arec, report$Brec, report$M)
    report <- c(report, ref_pt[1:3])

    Assessment@FMSY <- report$FMSY
    Assessment@MSY <- report$MSY
    Assessment@BMSY <- Assessment@SSBMSY <- Assessment@VBMSY <- report$BMSY
    Assessment@F_FMSY <- structure(report$F/report$FMSY, names = Year)
    Assessment@B_BMSY <- Assessment@SSB_SSBMSY <- Assessment@VB_VBMSY <- structure(report$B/report$BMSY, names = Yearplusone)
    Assessment@TMB_report <- report

    if(state_space) {
      Assessment@SE_Dev <- structure(as.list(SD, "Std. Error")$log_rec_dev, names = Year)
    }
    
    catch_eq <- function(Ftarget) {
      projection_cDD(Assessment, Ftarget = Ftarget, p_years = 1, p_sim = 1, obs_error = list(matrix(1, 1, 1), matrix(1, 1, 1)), 
                     process_error = matrix(1, 1, 1)) %>% slot("Catch") %>% as.vector()
    }
    Assessment@forecast <- list(per_recruit = ref_pt[[4]], catch_eq = catch_eq)
  }
  return(Assessment)
}

ref_pt_cDD <- function(TMB_data, Arec, Brec, M) {
  opt2 <- optimize(yield_fn_cDD, interval = c(0, 3), M = M, Kappa = TMB_data$Kappa, 
                   Winf = TMB_data$Winf, wk = TMB_data$wk, SR = TMB_data$SR_type, 
                   Arec = Arec, Brec = Brec)
  FMSY <- opt2$minimum
  MSY <- -1 * opt2$objective
  BMSY <- MSY/FMSY
  
  F_PR <- seq(0, 2.5 * FMSY, length.out = 100)
  yield <- lapply(F_PR, yield_fn_cDD, M = M, Kappa = TMB_data$Kappa, 
                  Winf = TMB_data$Winf, wk = TMB_data$wk, SR = TMB_data$SR_type, 
                  Arec = Arec, Brec = Brec, opt = FALSE)
  
  SPR <- vapply(yield, getElement, numeric(1), "SPR")
  YPR <- vapply(yield, getElement, numeric(1), "YPR")
  
  return(list(FMSY = FMSY, MSY = MSY, BMSY = BMSY,
              per_recruit = data.frame(FM = F_PR, SPR = SPR/SPR[1], YPR = YPR)))
}

yield_fn_cDD <- function(x, M, Kappa, Winf, wk, SR, Arec, Brec, opt = TRUE, log_trans = FALSE) {
  if(log_trans) {
    FMort <- exp(x)
  } else {
    FMort <- x
  } 
  Z <- FMort + M
  BPR <- (wk + Kappa * Winf/Z)/(Z+Kappa)
  if(SR == "BH") Req <- (Arec * BPR - 1)/Brec/BPR
  if(SR == "Ricker") Req <- log(Arec * BPR)/Brec/BPR
  Beq <- BPR * Req
  YPR <- FMort * BPR
  Yield <- FMort * Beq
  if(opt) {
    return(-1 * Yield)
  } else {
    return(c(SPR = BPR, Yield = Yield, YPR = YPR, B = Beq, R = Req))
  }
}


cDD_dynamic_SSB0 <- function(obj, par = obj$env$last.par.best, ...) {
  dots <- list(...)
  dots$data$C_hist <- rep(1e-8, dots$data$ny)
  par[names(par) == "F_equilibrium"] <- 0
  
  obj2 <- MakeADFun(data = dots$data, parameters = dots$params, map = dots$map, 
                    random = obj$env$random, DLL = "SAMtool", silent = TRUE)
  obj2$report(par)$B
}

