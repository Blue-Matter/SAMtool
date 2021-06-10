#' Delay - Difference Stock Assessment in TMB
#'
#' A simple delay-difference assessment model using a
#' time-series of catches and a relative abundance index and coded in TMB. The model
#' can be conditioned on either (1) effort and estimates predicted catch or (2) catch and estimates a predicted index.
#' In the state-space version \code{DD_SS}, recruitment deviations from the stock-recruit relationship are estimated.
#'
#' @param x An index for the objects in \code{Data} when running in closed loop simulation.
#' Otherwise, equals to 1 when running an assessment.
#' @param Data An object of class \linkS4class{Data}.
#' @param condition A string to indicate whether to condition the model on catch or effort (ratio of catch and index).
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
#' @param prior A named list (R0, h, M, and q) to provide the mean and standard deviations of prior distributions for those parameters. R0, index q, and M priors are
#' lognormal (provide the mean in normal space, SD in lognormal space). Beverton-Holt steepness uses a beta prior, while Ricker steepness uses a normal prior.
#' For index q, provide a matrix for nsurvey rows and 2 columns (for mean and SD), with NA in rows corresponding to indices without priors.
#' For all others, provide a length-2 vector for the mean and SD.
#' See vignette for full description.
#' @param fix_h Logical, whether to fix steepness to value in \code{Data@@steep} in the assessment model.
#' Automatically false if a prior is used.
#' @param fix_sd Logical, whether the standard deviation of the data in the likelihood (index for conditioning on catch or
#' catch for conditioning on effort). If \code{TRUE}, the SD is fixed to value provided in \code{start} (if provided), otherwise,
#'  value based on either \code{Data@@CV_Cat} or \code{Data@@CV_Ind}.
#' @param fix_tau Logical, the standard deviation of the recruitment deviations is fixed. If \code{TRUE},
#' tau is fixed to value provided in \code{start} (if provided), otherwise, equal to 1.
#' @param dep The initial depletion in the first year of the model. A tight prior is placed on the model objective function
#' to estimate the equilibrium exploitation rate that corresponds to the initial depletion. Due to this tight prior, this F
#' should not be considered to be an independent model parameter.
#' @param LWT A named list of likelihood weights. For \code{LWT$Index}, a vector of likelihood weights for each survey, while
#' for \code{LWT$MW} a numeric.
#' @param n_itF Integer, the number of iterations to solve F within an annual time step when conditioning on catch.
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
#' @param control A named list of parameters regarding optimization to be passed to
#' \code{\link[stats]{nlminb}}.
#' @param inner.control A named list of arguments for optimization of the random effects, which
#' is passed on to \code{\link[TMB]{newton}} via \code{\link[TMB]{MakeADFun}}.
#' @param ... Additional arguments (not currently used).
#' @return An object of \code{\linkS4class{Assessment}} containing objects and output from TMB.
#' @details 
#' For \code{start} (optional), a named list of starting values of estimates can be provided for:
#' \itemize{
#' \item \code{R0} Unfished recruitment. Otherwise, Data@@OM$R0[x] is used in closed-loop, and 400\% of mean catch otherwise.
#' \item \code{h} Steepness. Otherwise, Data@@steep[x] is used, or 0.9 if empty.
#' \item \code{M} Natural mortality. Otherwise, Data@@Mort[x] is used.
#' \item \code{k} Age of knife-edge maturity. By default, the age of 50\% maturity calculated from the slots in the Data object.
#' \item \code{Rho} Delay-difference rho parameter. Otherwise, calculated from biological parameters in the Data object.
#' \item \code{Alpha} Delay-difference alpha parameter. Otherwise, calculated from biological parameters in the Data object.
#' \item \code{q_effort} Scalar coefficient when conditioning on effort (to scale to F). Otherwise, 1 is the default.
#' \item \code{F_equilibrium} Equilibrium exploitation rate leading into first year of the model (to determine initial depletion). By default, 0.
#' \item \code{omega} Lognormal SD of the catch (observation error) when conditioning on effort. By default, Data@@CV_Cat[x].
#' \item \code{tau} Lognormal SD of the recruitment deviations (process error) for \code{DD_SS}. By default, Data@@sigmaR[x].
#' \item \code{sigma} Lognormal SD of the index (observation error) when conditioning on catch. By default, Data@@CV_Ind[x]. Not
#' used if multiple indices are used.
#' \item \code{sigma_W} Lognormal SD of the mean weight (observation error). By default, 0.1.
#' }
#' 
#' Multiple indices are supported in the model. Data@@Ind, Data@@VInd, and Data@@SpInd are all assumed to be biomass-based.
#' For Data@@AddInd, Data@@I_units are used to identify a biomass vs. abundance-based index.
#' 
#' Similar to many other assessment
#' models, the model depends on assumptions such as stationary productivity and
#' proportionality between the abundance index and real abundance.
#' Unsurprisingly the extent to which these assumptions are
#' violated tends to be the biggest driver of performance for this method.
#' @author T. Carruthers & Z. Siders. Zach Siders coded the TMB function.
#' @references
#' Carruthers, T, Walters, C.J,, and McAllister, M.K. 2012. Evaluating methods that classify
#' fisheries stock status using only fisheries catch data. Fisheries Research 119-120:66-79.
#'
#' Hilborn, R., and Walters, C., 1992. Quantitative Fisheries Stock Assessment: Choice,
#' Dynamics and Uncertainty. Chapman and Hall, New York.
#' @section Required Data:
#' \itemize{
#' \item \code{DD_TMB}: Cat, Ind, Mort, L50, vbK, vbLinf, vbt0, wla, wlb, MaxAge
#' \item \code{DD_SS}: Cat, Ind, Mort, L50, vbK, vbLinf, vbt0, wla, wlb, MaxAge
#' }
#' @section Optional Data:
#' \itemize{
#' \item \code{DD_TMB}: steep
#' \item \code{DD_SS}: steep, CV_Cat
#' }
#' @examples
#' \donttest{
#' #### Observation-error delay difference model
#' res <- DD_TMB(x = 3, Data = MSEtool::SimulatedData)
#'
#' # Provide starting values
#' start <- list(h = 0.95)
#' res <- DD_TMB(x = 3, Data = MSEtool::SimulatedData, start = start)
#'
#' summary(res@@SD) # Parameter estimates
#'
#' ### State-space version
#' ### Set recruitment variability SD = 0.3 (since fix_tau = TRUE)
#' res <- DD_SS(x = 3, Data = MSEtool::SimulatedData, start = list(tau = 0.3))
#' }
#' @seealso \link{plot.Assessment} \link{summary.Assessment} \link{retrospective} \link{profile} \link{make_MP}
#' @export
DD_TMB <- function(x = 1, Data, condition = c("catch", "effort"), AddInd = "B", SR = c("BH", "Ricker"), rescale = "mean1", MW = FALSE,
                   start = NULL, prior = list(), fix_h = TRUE, dep = 1, LWT = list(), nit_F = 3L, silent = TRUE, opt_hess = FALSE, n_restart = ifelse(opt_hess, 0, 1),
                   control = list(iter.max = 5e3, eval.max = 1e4), ...) {
  condition <- match.arg(condition)
  DD_(x = x, Data = Data, state_space = FALSE, condition = condition, AddInd = AddInd, SR = SR, rescale = rescale, 
      MW = MW, start = start, prior = prior, fix_h = fix_h, dep = dep, LWT = LWT, fix_sd = FALSE,
      fix_tau = TRUE, nit_F = nit_F, integrate = FALSE, silent = silent, opt_hess = opt_hess, n_restart = n_restart,
      control = control, inner.control = list(), ...)
}
class(DD_TMB) <- "Assess"


#' @rdname DD_TMB
#' @export
DD_SS <- function(x = 1, Data, condition = c("catch", "effort"), AddInd = "B", SR = c("BH", "Ricker"), rescale = "mean1", MW = FALSE,
                  start = NULL, prior = list(), fix_h = TRUE, fix_sd = FALSE, fix_tau = TRUE, dep = 1, LWT = list(),
                  nit_F = 3L, integrate = FALSE, silent = TRUE, opt_hess = FALSE, n_restart = ifelse(opt_hess, 0, 1),
                  control = list(iter.max = 5e3, eval.max = 1e4), inner.control = list(), ...) {
  condition <- match.arg(condition)
  DD_(x = x, Data = Data, state_space = TRUE, condition = condition, AddInd = AddInd, SR = SR, rescale = rescale, 
      MW = MW, start = start, prior = prior, fix_h = fix_h, dep = dep, LWT = LWT, fix_sd = fix_sd,
      fix_tau = fix_tau, nit_F = nit_F, integrate = integrate, silent = silent, opt_hess = opt_hess, n_restart = n_restart,
      control = control, inner.control = inner.control, ...)
}
class(DD_SS) <- "Assess"

#' @useDynLib SAMtool
DD_ <- function(x = 1, Data, state_space = FALSE, condition = c("catch", "effort"), AddInd = "B", SR = c("BH", "Ricker"), 
                rescale = "mean1", MW = FALSE, start = NULL, prior = list(), fix_h = TRUE, fix_sd = TRUE, fix_tau = TRUE, dep = 1, LWT = list(),
                nit_F = 3L, integrate = FALSE, silent = TRUE, opt_hess = FALSE, n_restart = ifelse(opt_hess, 0, 1),
                control = list(iter.max = 5e3, eval.max = 1e4), inner.control = list(), ...) {
  dependencies <- "Data@Cat, Data@Ind, Data@Mort, Data@L50, Data@vbK, Data@vbLinf, Data@vbt0, Data@wla, Data@wlb, Data@MaxAge"
  dots <- list(...)
  start <- lapply(start, eval, envir = environment())

  condition <- match.arg(condition)
  SR <- match.arg(SR)
  Winf <- Data@wla[x] * Data@vbLinf[x]^Data@wlb[x]
  age <- 1:Data@MaxAge
  la <- Data@vbLinf[x] * (1 - exp(-Data@vbK[x] * ((age - Data@vbt0[x]))))
  wa <- Data@wla[x] * la^Data@wlb[x]
  a50V <- iVB(Data@vbt0[x], Data@vbK[x], Data@vbLinf[x], Data@L50[x])
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
  Ind <- lapply(AddInd, Assess_I_hist, Data = Data, x = x, yind = yind)
  I_hist <- vapply(Ind, getElement, numeric(ny), "I_hist")
  I_sd <- vapply(Ind, getElement, numeric(ny), "I_sd")
  I_units <- vapply(Ind, getElement, numeric(1), "I_units")

  if(is.null(I_hist)) stop("No indices found.", call. = FALSE)
  nsurvey <- ncol(I_hist)

  if(condition == "effort") {
    if(nsurvey > 1) message("Only one index time series can be used when conditioning on effort.", call. = FALSE)
    E_hist <- C_hist/I_hist[, 1]
    if(any(is.na(E_hist))) stop("Missing values in catch and index in Data object.")
    E_rescale <- 1/mean(E_hist)
    E_hist <- E_hist * E_rescale
  } else {
    E_hist <- rep(1, length(yind))
  }
  
  if(MW) {
    if(!is.null(Data@Misc[[x]]$MW)) {
      MW_hist <- Data@Misc[[x]]$MW
    } else {
      MW_hist <- apply(Data@CAL[x, , ], 1, function(xx) {
        weighted.mean(Data@wla[x]*Data@CAL_mids^Data@wlb[x], xx, na.rm = TRUE)
      })
    }
    MW_hist[MW_hist <= 0] <- NA_real_
  } else {
    MW_hist <- rep(NA_real_, ny)
  }
  
  # Generate priors
  prior <- make_prior(prior, nsurvey, ifelse(SR == "BH", 1, 2), msg = FALSE)

  if(!is.null(start$k)) {
    k <- start$k
  } else {
    k <- ceiling(a50V)  # get age nearest to 50% vulnerability (ascending limb)
    k[k > Data@MaxAge/2] <- ceiling(Data@MaxAge/2)  # to stop stupidly high estimates of age at 50% vulnerability
  }
  if(!is.null(start$Rho)) {
    Rho <- start$Rho
  } else {
    Rho <- (wa[k + 2] - Winf)/(wa[k + 1] - Winf)
  }
  if(!is.null(start$Alpha)) {
    Alpha <- start$Alpha
  } else {
    Alpha <- Winf * (1 - Rho)
  }
  M <- Data@Mort[x]
  wk <- wa[k]

  if(rescale == "mean1") rescale <- 1/mean(C_hist)
  if(dep <= 0 || dep > 1) stop("Initial depletion (dep) must be between > 0 and <= 1.")
  if(!is.list(LWT)) {
    if(!is.null(LWT) && length(LWT) != nsurvey) stop("LWT needs to be a vector of length ", nsurvey)
    LWT <- list(Index = LWT)
    LWT$MW <- 1
  } else {
    if(is.null(LWT$Index)) LWT$Index <- rep(1, nsurvey)
    if(is.null(LWT$MW)) LWT$MW <- 1 
  }

  fix_sigma <- condition == "effort" || nsurvey > 1 || MW || fix_sd
  fix_omega <- condition == "catch" || MW || fix_sd
  data <- list(model = "DD", Alpha = Alpha, Rho = Rho, ny = ny, k = k,
               wk = wk, C_hist = C_hist, dep = dep, rescale = rescale, I_hist = I_hist, I_units = I_units, I_sd = I_sd,
               E_hist = E_hist, MW_hist = MW_hist, SR_type = SR, condition = condition, LWT = c(LWT$Index, LWT$MW),
               nsurvey = nsurvey, fix_sigma = as.integer(fix_sigma), nit_F = nit_F, state_space = as.integer(state_space),
               use_prior = prior$use_prior, prior_dist = prior$pr_matrix)
  LH <- list(LAA = la, WAA = wa, maxage = Data@MaxAge, A50 = k)

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
    if(!is.null(start$q_effort) && is.numeric(start$q_effort)) params$log_q_effort <- log(start$q_effort[1])
    if(!is.null(start$F_equilibrium) && is.numeric(start$F_equilibrium)) params$F_equilibrium <- start$F_equilibrium
    if(!is.null(start$omega) && is.numeric(start$omega)) params$log_omega <- log(start$omega[1])
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
  if(is.null(params$log_q_effort)) params$log_q_effort <- log(1)
  if(is.null(params$F_equilibrium)) params$F_equilibrium <- ifelse(dep < 1, 0.1, 0)
  if(is.null(params$log_omega)) {
    params$log_omega <- max(0.05, sdconv(1, Data@CV_Cat[x]), na.rm = TRUE) %>% log()
  }
  if(is.null(params[["log_sigma"]])) params$log_sigma <- max(0.05, sdconv(1, Data@CV_Ind[x]), na.rm = TRUE) %>% log()
  if(is.null(params[["log_sigma_W"]])) params$log_sigma_W <- log(0.1)
  if(is.null(params$log_tau)) {
    params$log_tau <- ifelse(is.na(Data@sigmaR[x]), 0.6, Data@sigmaR[x]) %>% log()
  }
  params$log_rec_dev = rep(0, ny)

  info <- list(Year = Year, data = data, params = params, I_hist = I_hist, LH = LH,
               rescale = rescale, control = control, inner.control = inner.control)
  if(condition == "effort") info$E_rescale <- E_rescale

  map <- list()
  if(condition == "catch") map$log_q_effort <- factor(NA)
  if(fix_h && !prior$use_prior[2]) map$transformed_h <- factor(NA)
  if(!prior$use_prior[3]) map$log_M <- factor(NA)
  if(dep == 1) map$F_equilibrium <- factor(NA)
  if(fix_omega) map$log_omega <- factor(NA)
  if(fix_sigma) map$log_sigma <- factor(NA)
  map$log_sigma_W <- factor(NA)
  if(fix_tau) map$log_tau <- factor(NA)
  if(!state_space) map$log_rec_dev <- factor(rep(NA, ny))

  random <- NULL
  if(integrate) random <- "log_rec_dev"
  
  obj <- MakeADFun(data = info$data, parameters = info$params, random = random,
                   map = map, hessian = TRUE, DLL = "SAMtool", inner.control = inner.control, silent = silent)

  mod <- optimize_TMB_model(obj, control, opt_hess, n_restart)
  opt <- mod[[1]]
  SD <- mod[[2]]
  report <- obj$report(obj$env$last.par.best)

  Yearplusone <- c(Year, max(Year) + 1)
  Yearplusk <- c(Year, max(Year) + 1:k)

  if(condition == "catch") {
    NLL_name <- paste0("Index_", 1:nsurvey)
  } else {
    NLL_name <- "Catch"
  }
  nll_report <- ifelse(is.character(opt), ifelse(integrate, NA, report$nll), opt$objective)
  report$dynamic_SSB0 <- DD_dynamic_SSB0(obj, data = info$data, params = info$params, map = map) %>% 
    structure(names = Yearplusone)
  Assessment <- new("Assessment", Model = ifelse(state_space, "DD_SS", "DD_TMB"),
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
                                    names = c("Total", NLL_name, "MW", "Dev", "Prior", "Penalty")),
                    info = info, obj = obj, opt = opt, SD = SD, TMB_report = report,
                    dependencies = dependencies)
  if(!MW) Assessment@NLL <- Assessment@NLL[names(Assessment@NLL) != "MW"]

  if(state_space) {
    Assessment@Dev <- structure(report$log_rec_dev, names = Year)
    Assessment@Dev_type <- "log-Recruitment deviations"
  }

  if(Assessment@conv) {
    ref_pt <- ref_pt_DD(info$data, report$Arec, report$Brec, report$M)
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
      projection_DD_TMB(Assessment, Ftarget = Ftarget, p_years = 1, p_sim = 1, obs_error = list(matrix(1, 1, 1), matrix(1, 1, 1)), 
                        process_error = matrix(1, 1, 1)) %>% slot("Catch") %>% as.vector()
    }
    Assessment@forecast <- list(per_recruit = ref_pt[[4]], catch_eq = catch_eq)
  }
  return(Assessment)
}


ref_pt_DD <- function(TMB_data, Arec, Brec, M) {
  opt2 <- optimize(yield_fn_DD, interval = c(1e-4, 3), M = M, Alpha = TMB_data$Alpha,
                   Rho = TMB_data$Rho, wk = TMB_data$wk, SR = TMB_data$SR_type, Arec = Arec, Brec = Brec)
  opt3 <- yield_fn_DD(opt2$minimum, M = M, Alpha = TMB_data$Alpha, opt = FALSE,
                      Rho = TMB_data$Rho, wk = TMB_data$wk, SR = TMB_data$SR_type, Arec = Arec, Brec = Brec)
  FMSY <- opt2$minimum
  MSY <- -1 * opt2$objective
  BMSY <- opt3["B"]
  
  F_PR <- seq(0, 2.5 * FMSY, length.out = 100)
  yield <- lapply(F_PR, 
                  yield_fn_DD, M = M, Alpha = TMB_data$Alpha,
                  Rho = TMB_data$Rho, wk = TMB_data$wk, SR = TMB_data$SR_type, 
                  Arec = Arec, Brec = Brec, opt = FALSE)
  SPR <- vapply(yield, getElement, numeric(1), "SPR")
  YPR <- vapply(yield, getElement, numeric(1), "YPR")
  
  return(list(FMSY = FMSY, MSY = MSY, BMSY = BMSY, 
              per_recruit = data.frame(FM = F_PR, SPR = SPR/SPR[1], YPR = YPR)))
}

yield_fn_DD <- function(x, M, Alpha, Rho, wk, SR, Arec, Brec, opt = TRUE, log_trans = FALSE) {
  if(log_trans) {
    FF <- exp(x)
  } else {
    FF <- x
  }
  S0 <- exp(-M)
  SS <- S0 * exp(-FF)
  Spr <- (SS * Alpha/(1 - SS) + wk)/(1 - Rho * SS)
  if(SR == "BH") Req <- (Arec * Spr - 1)/(Brec * Spr)
  if(SR == "Ricker") Req <- log(Arec * Spr)/(Brec * Spr)
  Beq <- Spr * Req
  Ypr <- FF * Spr * (1 - SS) / (FF + M)
  Yield <- Ypr * Req
  if(opt) {
    return(-1 * Yield)
  } else {
    return(c(SPR = Spr, Yield = Yield, YPR = Ypr, B = Beq, R = Req))
  }
}

DD_dynamic_SSB0 <- function(obj, par = obj$env$last.par.best, ...) {
  dots <- list(...)
  
  if(dots$data$condition == "catch") {
    dots$data$C_hist <- rep(1e-8, dots$data$ny)
  } else if(dots$data$condition == "effort") {
    par[names(par) == "log_q_effort"] <- -1e8
  }
  par[names(par) == "F_equilibrium"] <- 0
  
  obj2 <- MakeADFun(data = dots$data, parameters = dots$params, map = dots$map, 
                    random = obj$env$random, DLL = "SAMtool", silent = TRUE)
  obj2$report(par)$B
}

