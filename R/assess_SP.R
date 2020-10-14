#' Surplus production model with FMSY and MSY as leading parameters
#'
#' A surplus production model that uses only a time-series of catches and a relative abundance index
#' and coded in TMB. The base model, \code{SP}, is conditioned on catch and estimates a predicted index.
#' Continuous surplus production and fishing is modeled with sub-annual time steps which should approximate
#' the behavior of ASPIC (Prager 1994). The Fox model, \code{SP_Fox}, fixes BMSY/K = 0.37 (1/e).
#' The state-space version, \code{SP_SS} estimates annual deviates in biomass. An option allows for setting a
#' prior for the intrinsic rate of increase.
#' The function for the \code{spict} model (Pedersen and Berg, 2016) is available in \link[DLMtool]{DLMextra}.
#'
#' @param x An index for the objects in \code{Data} when running in \link[DLMtool]{runMSE}.
#' Otherwise, equals to 1 When running an assessment interactively.
#' @param Data An object of class Data.
#' @param AddInd A vector of integers or character strings indicating the indices to be used in the model. Integers assign the index to
#' the corresponding index in Data@@AddInd, "B" (or 0) represents total biomass in Data@@Ind, "VB" represents vulnerable biomass in
#' Data@@VInd, and "SSB" represents spawning stock biomass in Data@@SpInd.
#' @param rescale A multiplicative factor that rescales the catch in the assessment model, which
#' can improve convergence. By default, \code{"mean1"} scales the catch so that time series mean is 1, otherwise a numeric.
#' Output is re-converted back to original units.
#' @param start Optional list of starting values. Entries can be expressions that are evaluated in the function. See details.
#' @param fix_dep Logical, whether to fix the initial depletion (ratio of biomass to carrying capacity in the
#' first year of the model). If \code{TRUE}, uses the value in \code{start}, otherwise equal to 1
#' (unfished conditions).
#' @param fix_n Logical, whether to fix the exponent of the production function. If \code{TRUE},
#' uses the value in \code{start}, otherwise equal to \code{n = 2}, where the biomass at MSY
#' is half of carrying capacity.
#' @param fix_sigma Logical, whether the standard deviation of the index is fixed. If \code{TRUE},
#' sigma is fixed to value provided in \code{start} (if provided), otherwise, value based on \code{Data@@CV_Ind}.
#' @param fix_tau Logical, the standard deviation of the biomass deviations is fixed. If \code{TRUE},
#' tau is fixed to value provided in \code{start} (if provided), otherwise, equal to 0.1.
#' @param early_dev Character string describing the years for which biomass deviations are estimated in \code{SP_SS}.
#' By default, deviations are estimated in each year of the model (\code{"all"}), while deviations could also be estimated
#' once index data are available (\code{"index"}).
#' @param LWT A vector of likelihood weights for each survey.
#' @param n_seas Integer, the number of seasons in the model for calculating continuous surplus production.
#' @param n_itF Integer, the number of iterations to solve F conditional on the observed catch given multiple seasons within an annual time step.
#' Ignored if \code{n_seas} = 1.
#' @param integrate Logical, whether the likelihood of the model integrates over the likelihood
#' of the biomass deviations (thus, treating it as a state-space variable).
#' @param use_r_prior Logical, whether a prior for the intrinsic rate of increase will be used in the model. See details.
#' @param r_reps If \code{use_r_prior = TRUE}, the number of samples of natural mortality and steepness for calculating the
#' mean and standard deviation of the r prior. To override and directly provide the r-prior mean and standard deviation, use the start list, e.g.
#' \code{start = list(r_prior = c(0.1, 0.05))} (mean of 0.1 and s.d. of 0.05).
#' @param SR_type If \code{use_r_prior = TRUE}, the stock-recruit relationship used to calculate unfished recruits per spawner at the origin
#' of spwaning biomass approaches zero. Used for the r prior.
#' @param silent Logical, passed to \code{\link[TMB]{MakeADFun}}, whether TMB
#' will print trace information during optimization. Used for dignostics for model convergence.
#' @param opt_hess Logical, whether the hessian function will be passed to \code{\link[stats]{nlminb}} during optimization
#' (this generally reduces the number of iterations to convergence, but is memory and time intensive and does not guarantee an increase
#' in convergence rate). Ignored if \code{integrate = TRUE}.
#' @param n_restart The number of restarts (calls to \code{\link[stats]{nlminb}}) in the optimization procedure, so long as the model
#' hasn't converged. The optimization continues from the parameters from the previous (re)start.
#' @param control A named list of parameters regarding optimization to be passed to
#' \code{\link[stats]{nlminb}}.
#' @param inner.control A named list of arguments for optimization of the random effects, which
#' is passed on to \link[TMB]{newton} via \code{\link[TMB]{MakeADFun}}.
#' @param ... For \code{SP_Fox}, additional arguments to pass to \code{SP}.
#' @details
#' To provide starting values for the \code{SP}, a named list can be provided for \code{FMSY},
#' \code{MSY}, \code{dep}, and \code{n} via the start argument (see example).
#'
#' For \code{SP_SS}, a start value can also be provided for \code{sigma} and \code{tau}, the standard deviation
#' of the index and log-biomass deviates, respectively. Default for tau is 0.1. Deviations are estimated beginning in the year when index
#' data are available.
#'
#' If \code{use_r_prior = TRUE}, \code{SP} and \code{SP_SS} will use a prior for the intrinsic rate of increase in the objective function.
#' A vector of length two can be passed in the \code{start} list for the mean and standard deviation of the prior (see example). The normal
#' distribution is used.
#'
#' If no values are provided, a prior is created using the Euler-Lotka method (Equation 15a of McAllister et al. 2001).
#' The Euler-Lotka method is modified to multiply the left-hand side of equation 15a by the alpha parameter of the
#' stock-recruit relationship (Stanley et al. 2009). Natural mortality and steepness are sampled in order to generate
#' a prior distribution for r. See \code{vignette("Surplus_production")} for more details.
#' @return An object of \code{\linkS4class{Assessment}} containing objects and output from TMB.
#' @note The model uses the Fletcher (1978) formulation and is parameterized with FMSY and MSY as
#' leading parameters. The default conditions assume unfished conditions in the first year of the time series
#' and a symmetric production function (n = 2).
#'
#' Tip: to create the Fox model (Fox 1970), just fix n = 1. See example.
#' @author Q. Huynh
#' @references
#' Fletcher, R. I. 1978. On the restructuring of the Pella-Tomlinson system. Fishery Bulletin 76:515:521.
#'
#' Fox, W.W. 1970. An exponential surplus-yield model for optimizing exploited fish populations. Transactions of the American Fisheries Society 99:80-88.
#'
#' McAllister, M.K., Pikitch, E.K., and Babcock, E.A. 2001. Using demographic methods to construct Bayesian priors
#' for the intrinsic rate of increase in the Schaefer model and implications for stock rebuilding. Can. J. Fish.
#' Aquat. Sci. 58: 1871-1890.
#'
#' Pedersen, M. W. and Berg, C. W. 2017. A stochastic surplus production model in continuous time. Fish and Fisheries. 18:226-243.
#'
#' Pella, J. J. and Tomlinson, P. K. 1969. A generalized stock production model. Inter-Am. Trop. Tuna Comm., Bull. 13:419-496.
#'
#' Prager, M. H. 1994. A suite of extensions to a nonequilibrium surplus-production model. Fishery Bulletin 92:374-389.
#'
#' Stanley, R.D., M. McAllister, P. Starr and N. Olsen. 2009. Stock assessment for bocaccio (Sebastes
#' paucispinis) in British Columbia waters. DFO Can. Sci. Advis. Sec. Res. Doc. 2009/055. xiv + 200 p.
#' @section Required Data:
#' \itemize{
#' \item \code{SP}: Cat, Ind
#' \item \code{SP_SS}: Cat, Ind
#' }
#' @section Optional Data:
#' \code{SP_SS}: CV_Ind
#' @examples
#' data(swordfish)
#'
#' #### Observation-error surplus production model
#' res <- SP(Data = swordfish)
#'
#' # Provide starting values, assume B/K = 0.875 in first year of model
#' # and symmetrical production curve (n = 2)
#' start <- list(dep = 0.875, n = 2)
#' res <- SP(Data = swordfish, start = start)
#'
#' \dontrun{
#' plot(res)
#' }
#' \donttest{
#' profile(res, FMSY = seq(0.1, 0.4, 0.01))
#' retrospective(res)
#' }
#'
#' #### State-space version
#' res_SS <- SP_SS(Data = swordfish, start = list(dep = 0.875, sigma = 0.1, tau = 0.1))
#'
#' \dontrun{
#' plot(res_SS)
#' }
#'
#' #### Fox model
#' res_Fox <- SP(Data = swordfish, start = list(n = 1), fix_n = TRUE)
#' res_Fox2 <- SP_Fox(Data = swordfish)
#'
#' #### SP with r_prior
#' res_prior <- SP(Data = SimulatedData, use_r_prior = TRUE)
#'
#' #### Pass an r_prior to the model with mean = 0.35, sd = 0.10
#' res_prior2 <- SP(Data = SimulatedData, use_r_prior = TRUE, start = list(r_prior = c(0.35, 0.10)))
#' @seealso \link{SP_production} \link{plot.Assessment} \link{summary.Assessment} \link{retrospective} \link{profile} \link{make_MP}
#' @export
SP <- function(x = 1, Data, AddInd = "B", rescale = "mean1", start = NULL, fix_dep = TRUE, fix_n = TRUE, LWT = NULL,
               n_seas = 4L, n_itF = 3L, use_r_prior = FALSE, r_reps = 1e2, SR_type = c("BH", "Ricker"),
               silent = TRUE, opt_hess = FALSE, n_restart = ifelse(opt_hess, 0, 1),
               control = list(iter.max = 5e3, eval.max = 1e4), ...) {
  SP_(x = x, Data = Data, AddInd = AddInd, state_space = FALSE, rescale = rescale, start = start, fix_dep = fix_dep, fix_n = fix_n, fix_sigma = TRUE,
      fix_tau = TRUE, LWT = LWT, n_seas = n_seas, n_itF = n_itF, use_r_prior = use_r_prior, r_reps = r_reps, SR_type = SR_type, integrate = FALSE,
      silent = silent, opt_hess = opt_hess, n_restart = n_restart, control = control, inner.control = list(), ...)
}
class(SP) <- "Assess"


#' @rdname SP
#' @export
SP_SS <- function(x = 1, Data, AddInd = "B", rescale = "mean1", start = NULL, fix_dep = TRUE, fix_n = TRUE, fix_sigma = TRUE,
                  fix_tau = TRUE, LWT = NULL, early_dev = c("all", "index"), n_seas = 4L, n_itF = 3L,
                  use_r_prior = FALSE, r_reps = 1e2, SR_type = c("BH", "Ricker"), integrate = FALSE,
                  silent = TRUE, opt_hess = FALSE, n_restart = ifelse(opt_hess, 0, 1),
                  control = list(iter.max = 5e3, eval.max = 1e4), inner.control = list(), ...) {
  SP_(x = x, Data = Data, AddInd = AddInd, state_space = TRUE, rescale = rescale, start = start, fix_dep = fix_dep, fix_n = fix_n, fix_sigma = fix_sigma,
      fix_tau = fix_tau, early_dev = early_dev, LWT = LWT, n_seas = n_seas, n_itF = n_itF, use_r_prior = use_r_prior, r_reps = r_reps,
      SR_type = SR_type, integrate = integrate, silent = silent, opt_hess = opt_hess, n_restart = n_restart,
      control = control, inner.control = inner.control, ...)
}
class(SP_SS) <- "Assess"

#' @rdname SP
#' @export
SP_Fox <- function(x = 1, Data, ...) {
  SP_args <- c(x = x, Data = Data, list(...))
  SP_args$start$n <- 1
  SP_args$fix_n <- TRUE

  do.call(SP, SP_args)
}
class(SP_Fox) <- "Assess"


#' @importFrom TMB MakeADFun
#' @importFrom stats nlminb
#' @useDynLib MSEtool
SP_ <- function(x = 1, Data, AddInd = "B", state_space = FALSE, rescale = "mean1", start = NULL, fix_dep = TRUE, fix_n = TRUE, fix_sigma = TRUE,
                fix_tau = TRUE, early_dev = c("all", "index"), LWT = NULL, n_seas = 4L, n_itF = 3L,
                use_r_prior = FALSE, r_reps = 1e2, SR_type = c("BH", "Ricker"), integrate = FALSE,
                silent = TRUE, opt_hess = FALSE, n_restart = ifelse(opt_hess, 0, 1),
                control = list(iter.max = 5e3, eval.max = 1e4), inner.control = list(), ...) {

  dependencies = "Data@Cat, Data@Ind"
  dots <- list(...)
  start <- lapply(start, eval, envir = environment())

  early_dev <- match.arg(early_dev)
  if(any(names(dots) == "yind")) {
    yind <- eval(dots$yind)
  } else {
    ystart <- which(!is.na(Data@Cat[x, ]))[1]
    yind <- ystart:length(Data@Cat[x, ])
  }
  Year <- Data@Year[yind]
  C_hist <- Data@Cat[x, yind]
  if(any(is.na(C_hist))) stop('Model is conditioned on complete catch time series, but there is missing catch.')
  ny <- length(C_hist)
  if(rescale == "mean1") rescale <- 1/mean(C_hist)

  Ind <- lapply(AddInd, Assess_I_hist, Data = Data, x = x, yind = yind)
  I_hist <- do.call(cbind, lapply(Ind, getElement, "I_hist"))
  I_sd <- do.call(cbind, lapply(Ind, getElement, "I_sd"))
  if(is.null(I_hist)) stop("No indices found.")
  nsurvey <- ncol(I_hist)

  if(state_space) {
    if(early_dev == "all") est_B_dev <- rep(1, ny)
    if(early_dev == "index") {
      first_year_index <- which(apply(I_hist, 1, function(x) any(!is.na(x))))[1]
      est_B_dev <- ifelse(1:ny < first_year_index, 0, 1)
    }
  } else {
    if(nsurvey == 1 && (AddInd == 0 | AddInd == "B")) {
      fix_sigma <- FALSE # Override: estimate sigma if there's a single survey
    }
    est_B_dev <- rep(0, ny)
  }

  if(is.null(LWT)) LWT <- rep(1, nsurvey)
  if(length(LWT) != nsurvey) stop("LWT needs to be a vector of length ", nsurvey)
  data <- list(model = "SP", C_hist = C_hist, rescale = rescale, I_hist = I_hist, I_sd = I_sd, I_lambda = LWT,
               fix_sigma = as.integer(fix_sigma), nsurvey = nsurvey, ny = ny,
               est_B_dev = est_B_dev, nstep = n_seas, dt = 1/n_seas, nitF = n_itF)

  if(use_r_prior) {
    if(!is.null(start$r_prior) && length(start$r_prior) == 2) {
      rp <- data$r_prior <- start$r_prior
    } else {
      rp <- r_prior_fn(x, Data, r_reps = r_reps, SR_type = SR_type)
      data$r_prior <- c(mean(rp), max(sd(rp), 0.1 * mean(rp)))
    }
  } else {
    rp <- data$r_prior <- c(0, 0)
  }

  params <- list()
  if(!is.null(start)) {
    if(!is.null(start$FMSY) && is.numeric(start$FMSY)) params$log_FMSY <- log(start$FMSY[1])
    if(!is.null(start$MSY) && is.numeric(start$MSY)) params$MSYx <- log(start$MSY[1])
    if(!is.null(start$dep) && is.numeric(start$dep)) params$log_dep <- log(start$dep[1])
    if(!is.null(start$n) && is.numeric(start$n)) params$log_n <- log(start$n[1])
    if(!is.null(start$sigma) && is.numeric(start$sigma)) params$log_sigma <- log(start$sigma)
    if(!is.null(start$tau) && is.numeric(start$tau)) params$log_tau <- log(start$tau[1])
  }
  if(is.null(params$log_FMSY)) params$log_FMSY <- ifelse(is.na(Data@Mort[x]), 0.2, 0.5 * Data@Mort[x]) %>% log()
  if(is.null(params$MSYx)) params$MSYx <- mean(3 * C_hist * rescale) %>% log()
  if(is.null(params$log_dep)) params$log_dep <- log(1)
  if(is.null(params$log_n)) params$log_n <- log(2)
  if(is.null(params$log_sigma)) params$log_sigma <- rep(log(0.05), nsurvey)
  if(is.null(params$log_tau)) params$log_tau <- log(0.1)
  params$log_B_dev <- rep(0, ny)

  map <- list()
  if(fix_dep) map$log_dep <- factor(NA)
  if(fix_n) map$log_n <- factor(NA)
  if(fix_sigma) map$log_sigma <- factor(rep(NA, nsurvey))
  if(fix_tau) map$log_tau <- factor(NA)
  if(any(!est_B_dev)) map$log_B_dev <- factor(ifelse(est_B_dev, 1:sum(est_B_dev), NA))

  random <- NULL
  if(integrate) random <- "log_B_dev"

  info <- list(Year = Year, data = data, params = params, rp = rp, control = control, inner.control = inner.control)

  obj <- MakeADFun(data = info$data, parameters = info$params, hessian = TRUE,
                   map = map, random = random, DLL = "MSEtool", silent = silent)
  mod <- optimize_TMB_model(obj, control, opt_hess, n_restart)
  opt <- mod[[1]]
  SD <- mod[[2]]
  report <- obj$report(obj$env$last.par.best)

  Yearplusone <- c(Year, max(Year) + 1)

  nll_report <- ifelse(is.character(opt), ifelse(integrate, NA, report$nll), opt$objective)
  Assessment <- new("Assessment", Model = ifelse(state_space, "SP_SS", "SP"), Name = Data@Name, conv = !is.character(SD) && SD$pdHess,
                    FMSY = report$FMSY, MSY = report$MSY, BMSY = report$BMSY, VBMSY = report$BMSY,
                    B0 = report$K, VB0 = report$K, FMort = structure(report$F, names = Year),
                    F_FMSY = structure(report$F/report$FMSY, names = Year),
                    B = structure(report$B, names = Yearplusone),
                    B_BMSY = structure(report$B/report$BMSY, names = Yearplusone),
                    B_B0 = structure(report$B/report$K, names = Yearplusone),
                    VB = structure(report$B, names = Yearplusone),
                    VB_VBMSY = structure(report$B/report$BMSY, names = Yearplusone),
                    VB_VB0 = structure(report$B/report$K, names = Yearplusone),
                    SSB = structure(report$B, names = Yearplusone),
                    SSB_SSBMSY = structure(report$B/report$BMSY, names = Yearplusone),
                    SSB_SSB0 = structure(report$B/report$K, names = Yearplusone),
                    Obs_Catch = structure(C_hist, names = Year), Obs_Index = structure(I_hist, dimnames = list(Year, paste0("Index_", 1:nsurvey))),
                    Catch = structure(report$Cpred, names = Year), Index = structure(report$Ipred, dimnames = list(Year, paste0("Index_", 1:nsurvey))),
                    NLL = structure(c(nll_report, report$nll_comp, report$penalty, report$prior),
                                    names = c("Total", paste0("Index_", 1:nsurvey), "Dev", "Penalty", "Prior")),
                    info = info, obj = obj, opt = opt, SD = SD, TMB_report = report,
                    dependencies = dependencies)

  if(state_space) {
    Assessment@Dev <- structure(report$log_B_dev, names = Year)
    Assessment@Dev_type <- "log-Biomass deviations"
    Assessment@NLL <- structure(c(nll_report, report$nll_comp, report$penalty, report$prior),
                                names = c("Total", paste0("Index_", 1:nsurvey), "Dev", "Penalty", "Prior"))
  } else {
    Assessment@NLL <- structure(c(nll_report, report$nll_comp[1:nsurvey], report$penalty, report$prior),
                                names = c("Total", paste0("Index_", 1:nsurvey), "Penalty", "Prior"))
  }

  if(Assessment@conv) {
    if(state_space) {
      if(integrate) {
        SE_Dev <- ifelse(est_B_dev, sqrt(SD$diag.cov.random), 0)
      } else {
        SE_Dev <- ifelse(est_B_dev, sqrt(diag(SD$cov.fixed)[names(SD$par.fixed) =="log_B_dev"]), 0)
      }
      Assessment@SE_Dev <- structure(SE_Dev, names = Year)
    }
    Assessment@SE_FMSY <- SD$sd[names(SD$value) == "FMSY"]
    Assessment@SE_MSY <- SD$sd[names(SD$value) == "MSY"]
    Assessment@SE_F_FMSY_final <- SD$sd[names(SD$value) == "F_FMSY_final"]
    Assessment@SE_B_BMSY_final <- SD$sd[names(SD$value) == "B_BMSY_final"]
    Assessment@SE_B_B0_final <- SD$sd[names(SD$value) == "B_K_final"]
    Assessment@SE_VB_VBMSY_final <- SD$sd[names(SD$value) == "B_BMSY_final"]
    Assessment@SE_VB_VB0_final <- SD$sd[names(SD$value) == "B_K_final"]
  }
  return(Assessment)
}


r_prior_fn <- function(x = 1, Data, r_reps = 1e2, SR_type = c("BH", "Ricker"), seed = x) {
  SR_type <- match.arg(SR_type)

  set.seed(x)
  M <- trlnorm(r_reps, Data@Mort[x], Data@CV_Mort[x])
  steep <- sample_steepness3(r_reps, Data@steep[x], Data@CV_steep[x], SR_type)

  max_age <- Data@MaxAge
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

  log_r <- vapply(1:r_reps, function(y) uniroot(Euler_Lotka_fn, c(-6, 2), M = M[y], h = steep[y], weight = Wa,
                                                mat = mat_age, maxage = max_age, SR_type = SR_type)$root, numeric(1))
  return(exp(log_r))
}

Euler_Lotka_fn <- function(log_r, M, h, weight, mat, maxage, SR_type) {
  M <- rep(M, maxage)
  surv0 <- exp(-M)
  NPR <- c(1, cumprod(surv0[1:maxage-1]))
  NPR[maxage] <- NPR[maxage]/(1 - surv0[maxage])

  SBPR <- sum(NPR * weight * mat)
  CR <- ifelse(SR_type == "BH", 4*h/(1-h), (5*h)^1.25)
  R_per_S <- CR/SBPR

  EL <- R_per_S * sum(NPR * weight * mat * exp(-exp(log_r) * c(1:maxage)))
  return(EL - 1)
}
