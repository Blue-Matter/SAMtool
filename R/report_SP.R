summary_SP <- function(Assessment, state_space = FALSE) {
  assign_Assessment_slots(Assessment)

  current_status <- data.frame(Value = c(F_FMSY[length(F_FMSY)], B_BMSY[length(B_BMSY)],
                                         B_B0[length(B_B0)]))
  rownames(current_status) <- c("F/FMSY", "B/BMSY", "B/B0")

  Value <- numeric(0)
  Description <- character(0)
  rownam <- character(0)
  if("log_dep" %in% names(obj$env$map)) {
    Value <- c(Value, TMB_report$dep)
    Description <- c(Description, "Initial depletion")
    rownam <- c(rownam, "dep")
  }
  if("log_n" %in% names(obj$env$map)) {
    Value <- c(Value, TMB_report$n)
    Description <- c(Description, "Production exponent")
    rownam <- c(rownam, "n")
  }
  if(state_space && "log_tau" %in% names(obj$env$map)) {
    Value <- c(Value, TMB_report$tau)
    Description <- c(Description, "Biomass deviation SD (log-space)")
    rownam <- c(rownam, "tau")
  }
  if(length(Value) == 0) input_parameters <- data.frame() else {
    input_parameters <- data.frame(Value = Value, Description = Description, stringsAsFactors = FALSE)
    rownames(input_parameters) <- rownam
  }

  derived <- data.frame(Value = c(TMB_report$r, TMB_report$K, TMB_report$BMSY, TMB_report$BMSY/TMB_report$K),
                        Description = c("Intrinsic rate of population increase", "Carrying capacity",
                                        "Biomass at MSY", "Depletion at MSY"),
                        stringsAsFactors = FALSE)
  rownames(derived) <- c("r", "K", "BMSY", "BMSY/B0")

  model_estimates <- sdreport_int(SD)
  if(!is.character(model_estimates)) {
    rownames(model_estimates)[rownames(model_estimates) == "log_B_dev"] <- paste0("log_B_dev_", names(FMort)[as.logical(obj$env$data$est_B_dev)])
  }

  model_name <- "Surplus Production"
  if(state_space) model_name <- paste(model_name, "(State-Space)")
  output <- list(model = model_name, current_status = current_status,
                 input_parameters = input_parameters, derived_quantities = derived,
                 model_estimates = model_estimates,
                 log_likelihood = matrix(NLL, ncol = 1, dimnames = list(names(NLL), "Neg.LL")))
  return(output)
}

rmd_SP <- function(Assessment, state_space = FALSE, ...) {
  if(state_space) {
    ss <- rmd_summary("Surplus Production (State-Space)")
  } else ss <- rmd_summary("Surplus Production")

  # Data section
  data_section <- c(rmd_data_timeseries("Catch", header = "## Data\n"),
                    rmd_data_timeseries("Index", is_matrix = is.matrix(Assessment@Obs_Index), nsets = ncol(Assessment@Obs_Index)))

  # Assessment
  #### Pars and Fit
  assess_fit <- c(rmd_FMSY(header = "## Assessment {.tabset}\n### Estimates and Model Fit\n"), rmd_MSY(),
                  rmd_F_FMSY_terminal(), rmd_B_BMSY_terminal(), rmd_B_B0_terminal(),
                  rmd_assess_fit_series(nsets = ncol(Assessment@Index)),
                  rmd_assess_fit("Catch", "catch", match = TRUE))

  if(state_space) {
    assess_fit2 <- c(rmd_residual("Dev", fig.cap = "Time series of biomass deviations.", label = Assessment@Dev_type),
                     rmd_residual("Dev", "SE_Dev", fig.cap = "Time series of biomass deviations with 95% confidence intervals.",
                                  label = Assessment@Dev_type, conv_check = TRUE))
    assess_fit <- c(assess_fit, assess_fit2)
  }

  #### Time Series
  ts_output <- c(rmd_F(header = "### Time Series Output\n"), rmd_F_FMSY(FALSE), rmd_B(), rmd_B_BMSY(FALSE),
                 rmd_B_B0(FALSE), rmd_Kobe("B_BMSY", xlab = "expression(B/B[MSY])", conv_check = FALSE))

  productivity <- c(rmd_yield_F("SP", FALSE, header = "### Productivity\n"), rmd_yield_depletion("SP", FALSE), rmd_sp(FALSE))

  return(c(ss, data_section, assess_fit, ts_output, productivity))
}



profile_likelihood_SP <- function(Assessment, ...) {
  dots <- list(...)
  if(!"FMSY" %in% names(dots) && !"MSY" %in% names(dots)) stop("Sequence of neither FMSY nor MSY was found. See help file.")
  if(!is.null(dots$FMSY)) FMSY <- dots$FMSY else {
    FMSY <- Assessment@FMSY
    profile_par <- "MSY"
  }
  if(!is.null(dots$MSY)) MSY <- dots$MSY else {
    MSY <- Assessment@MSY
    profile_par <- "FMSY"
  }

  map <- Assessment@obj$env$map
  params <- Assessment@info$params

  profile_grid <- expand.grid(FMSY = FMSY, MSY = MSY)
  joint_profile <- !exists("profile_par")

  profile_fn <- function(i, Assessment, params, map) {
    params$log_FMSY <- log(profile_grid[i, 1])
    params$MSYx <- log(profile_grid[i, 2] * Assessment@obj$env$data$rescale)

    if(joint_profile && length(Assessment@obj$par) == 2) {
      nll <- Assessment@obj$fn(x = c(params$log_FMSY, params$MSYx))
    } else {
      if(joint_profile) map$MSYx <- map$log_FMSY <- factor(NA) else {
        if(profile_par == "MSY") map$MSYx <- factor(NA) else map$log_FMSY <- factor(NA)
      }
      obj2 <- MakeADFun(data = Assessment@info$data, parameters = params, map = map,
                        random = Assessment@obj$env$random, DLL = "MSEtool", silent = TRUE)
      opt2 <- optimize_TMB_model(obj2, Assessment@info$control)[[1]]
      if(!is.character(opt2)) nll <- opt2$objective else nll <- NA
    }
    return(nll)
  }
  nll <- vapply(1:nrow(profile_grid), profile_fn, numeric(1), Assessment = Assessment, params = params, map = map) - Assessment@opt$objective
  profile_grid$nll <- nll

  if(joint_profile) {
    pars <- c("FMSY", "MSY")
    MLE <- vapply(pars, function(x, y) slot(y, x), y = Assessment, numeric(1))
  } else {
    pars <- profile_par
    MLE <- slot(Assessment, pars)
  }

  output <- new("prof", Model = Assessment@Model, Name = Assessment@Name, Par = pars, MLE = MLE, grid = profile_grid)
  return(output)
}




retrospective_SP <- function(Assessment, nyr, state_space = FALSE) {
  assign_Assessment_slots(Assessment)
  ny <- info$data$ny
  Year <- info$Year
  Year <- c(Year, max(Year) + 1)

  # Array dimension: Retroyr, Year, ts
  # ts includes: F, F/FMSY, B, B/BMSY, B/B0
  retro_ts <- array(NA, dim = c(nyr + 1, ny + 1, 5))
  TS_var <- c("F", "F_FMSY", "B", "B_BMSY", "B_B0")
  dimnames(retro_ts) <- list(Peel = 0:nyr, Year = Year, Var = TS_var)

  retro_est <- array(NA, dim = c(nyr + 1, length(SD$par.fixed[names(SD$par.fixed) != "log_B_dev"]), 2))
  dimnames(retro_est) <- list(Peel = 0:nyr, Var = names(SD$par.fixed)[names(SD$par.fixed) != "log_B_dev"],
                              Value = c("Estimate", "Std. Error"))

  lapply_fn <- function(i, info, obj, state_space) {
    ny_ret <- ny - i
    info$data$ny <- ny_ret
    info$data$C_hist <- info$data$C_hist[1:ny_ret]
    info$data$I_hist <- info$data$I_hist[1:ny_ret, , drop = FALSE]
    info$data$I_sd <- info$data$I_sd[1:ny_ret, , drop = FALSE]

    map <- obj$env$map
    if(state_space) {
      info$data$est_B_dev <- info$data$est_B_dev[1:ny_ret]
      info$params$log_B_dev <- rep(0, ny_ret)
      map$log_B_dev <- obj$env$map$log_B_dev[1:ny_ret]
    }

    obj2 <- MakeADFun(data = info$data, parameters = info$params, map = map, random = obj$env$random,
                      DLL = "MSEtool", silent = TRUE)
    mod <- optimize_TMB_model(obj2, info$control)
    opt2 <- mod[[1]]
    SD <- mod[[2]]

    if(!is.character(opt2) && !is.character(SD)) {
      report <- obj2$report(obj2$env$last.par.best)

      FMort <- c(report$F, rep(NA, 1 + i))
      F_FMSY <- FMort/report$FMSY
      B <- c(report$B, rep(NA, i))
      B_BMSY <- B/report$BMSY
      B_B0 <- B/report$K

      retro_ts[i+1, , ] <<- cbind(FMort, F_FMSY, B, B_BMSY, B_B0)

      sumry <- summary(SD, "fixed")
      sumry <- sumry[rownames(sumry) != "log_B_dev", drop = FALSE]
      retro_est[i+1, , ] <<- sumry

      return(SD$pdHess)
    }
    return(FALSE)
  }

  conv <- vapply(0:nyr, lapply_fn, logical(1), info = info, obj = obj, state_space = state_space)
  if(any(!conv)) warning("Peels that did not converge: ", paste0(which(!conv) - 1, collapse = " "))

  retro <- new("retro", Model = Assessment@Model, Name = Assessment@Name, TS_var = TS_var, TS = retro_ts,
               Est_var = dimnames(retro_est)[[2]], Est = retro_est)
  attr(retro, "TS_lab") <- c("Fishing mortality", expression(F/F[MSY]), "Biomass", expression(B/B[MSY]), expression(B/B[0]))

  return(retro)
}


summary_SP_SS <- function(Assessment) summary_SP(Assessment, TRUE)

rmd_SP_SS <- function(Assessment, ...) rmd_SP(Assessment, TRUE, ...)

profile_likelihood_SP_SS <- profile_likelihood_SP

retrospective_SP_SS <- function(Assessment, nyr) retrospective_SP(Assessment, nyr, TRUE)





plot_yield_SP <- function(data = NULL, report, fmsy, msy, xaxis = c("F", "Biomass", "Depletion"), relative_yaxis = FALSE) {
  BKratio <- seq(0, 1, 0.01)

  K <- report$K
  n <- report$n
  BMSY <- report$BMSY

  if(n == 1) {
    Yield <- ifelse(BKratio == 0, 0, -exp(1) * msy * BKratio * log(BKratio))
  } else {
    gamma.par <- n^(n/(n-1))/(n-1)
    Yield <- gamma.par * msy * (BKratio - BKratio^n)
  }

  Biomass <- BKratio * K
  F.vector <- Yield/Biomass

  if(relative_yaxis) {
    Yield <- Yield/max(Yield)
    ylab <- "Relative Equilibrium Yield"
  } else ylab <- "Equilibrium Yield"

  if(xaxis == "F") {
    plot(F.vector, Yield, typ = 'l', xlab = "Fishing Mortality F", ylab = ylab)
    segments(x0 = fmsy, y0 = 0, y1 = max(Yield), lty = 2)
    segments(x0 = 0, y0 = max(Yield), x1 = fmsy, lty = 2)
    abline(h = 0, col = 'grey')
  }

  if(xaxis == "Biomass") {
    plot(Biomass, Yield, typ = 'l', xlab = "Biomass", ylab = ylab)
    segments(x0 = BMSY, y0 = 0, y1 = max(Yield), lty = 2)
    segments(x0 = 0, y0 = max(Yield), x1 = BMSY, lty = 2)
    abline(h = 0, col = 'grey')
  }

  if(xaxis == "Depletion") {
    plot(BKratio, Yield, typ = 'l', xlab = expression(B/B[0]), ylab = ylab)
    segments(x0 = BMSY/K, y0 = 0, y1 = max(Yield), lty = 2)
    segments(x0 = 0, y0 = max(Yield), x1 = BMSY/K, lty = 2)
    abline(h = 0, col = 'grey')
  }
  invisible(data.frame(F = F.vector, Yield = Yield, B = Biomass, B_B0 = BKratio))
}

#' Find the production parameter based on depletion that produces MSY
#'
#' For surplus production models, this function returns the production exponent n corresponding
#' to BMSY/K (Fletcher 1978).
#'
#' @param depletion The hypothesized depletion that produces MSY.
#' @param figure Local, plots figure of production function as a function of depletion (B/K)
#'
#' @author Q. Huynh
#' @references
#' Fletcher, R. I. 1978. On the restructuring of the Pella-Tomlinson system. Fishery Bulletin 76:515:521.
#' @note May be useful for parameterizing \code{n} in \link{SP} and \link{SP_SS}.
#' @examples SP_production(0.5)
#' @return The production function exponent n (numeric).
#' @importFrom stats uniroot
#' @examples
#' SP_production(0.5)
#' @seealso \link{SP} \link{SP_SS}
#' @export SP_production
SP_production <- function(depletion, figure = TRUE) {

  if(length(depletion) > 1) {
    depletion <- depletion[1]
    message(paste("Function is not vectorized. Depletion value of", depletion, "is used."))
  }
  if(depletion <= 0 || depletion >= 1) stop(paste("Proposed depletion =", depletion, "but value must be between 0 and 1."))

  calc_depletion <- function(n) {
    depletion_MSY <- if(n==1) 1/exp(1) else n^(1/(1-n))
    return(depletion_MSY)
  }
  n_solver <- function(x) calc_depletion(x) - depletion
  get_n <- uniroot(f = n_solver, interval = c(0, 1e3))
  n_answer <- round(get_n$root, 3)

  if(figure) {
    fmsy <- 0.1
    msy <- fmsy * depletion
    plot_yield_SP(report = list(n = n_answer, BMSY = depletion, K = 1), fmsy = fmsy,
                  msy = msy, xaxis = "Depletion", relative_yaxis = TRUE)
    title(paste0("Production exponent n = ", n_answer))
  }
  return(n_answer)
}














