
summary_cDD <- function(Assessment, state_space = FALSE) {
  assign_Assessment_slots(Assessment)

  if (conv) current_status <- c(F_FMSY[length(F_FMSY)], B_BMSY[length(B_BMSY)], B_B0[length(B_B0)])
  else current_status <- c(NA, NA, B_B0[length(B_B0)])
  current_status <- data.frame(Value = current_status)
  rownames(current_status) <- c("F/FMSY", "B/BMSY", "B/B0")

  Value <- as.numeric(c(info$data$Kappa, info$data$k, info$data$wk, info$data$Winf))
  Description <- c("Weight coefficient", "Age of knife-edge selectivity", 
                   "Weight at age k", "Asymptotic weight")
  rownam <- c("Kappa", "k", "w_k", "Winf")
  if (state_space && "log_tau" %in% names(obj$env$map)) {
    Value <- c(Value, TMB_report$tau)
    Description <- c(Description, "log-Recruitment deviation SD")
    rownam <- c(rownam, "tau")
  }
  if ("transformed_h" %in% names(obj$env$map)) {
    Value <- c(Value, h)
    Description <- c(Description, "Stock-recruit steepness")
    rownam <- c(rownam, "h")
  }
  if ("log_M" %in% names(obj$env$map)) {
    Value <- c(Value, TMB_report$M)
    Description <- c(Description, "Natural mortality")
    rownam <- c(rownam, "M")
  }
  if (any(info$data$MW_hist > 0, na.rm = TRUE) && "log_sigma_W" %in% names(obj$env$map)) {
    Value <- c(Value, TMB_report$sigma_W)
    Description <- c(Description, "Mean weight SD")
    rownam <- c(rownam, "sigma_W")
  }
  input_parameters <- data.frame(Value = Value, Description = Description, stringsAsFactors = FALSE)
  rownames(input_parameters) <- rownam

  if (conv) derived <- c(B0, N0, MSY, FMSY, BMSY, BMSY/B0)
  else derived <- rep(NA, 6)
  derived <- data.frame(Value = derived,
                        Description = c("Unfished biomass", "Unfished abundance", "Maximum sustainable yield (MSY)",
                                        "Fishing mortality at MSY", "Biomass at MSY", "Depletion at MSY"),
                        stringsAsFactors = FALSE)
  rownames(derived) <- c("B0", "N0", "MSY", "FMSY", "BMSY", "BMSY/B0")

  model_estimates <- sdreport_int(SD)
  if (!is.character(model_estimates)) {
    rownames(model_estimates)[rownames(model_estimates) == "log_rec_dev"] <- paste0("log_rec_dev_", names(Dev))
  }

  model_name <- "Continuous Delay-Differential"
  if (state_space) model_name <- paste(model_name, "(State-Space)")
  output <- list(model = model_name, current_status = current_status, input_parameters = input_parameters,
                 derived_quantities = derived, model_estimates = model_estimates,
                 log_likelihood = matrix(NLL, ncol = 1, dimnames = list(names(NLL), "Neg.LL")))
  return(output)
}



rmd_cDD <- function(Assessment, state_space = FALSE, ...) {
  if (state_space) {
    ss <- rmd_summary("Continuous Delay-Differential (State-Space)")
  } else ss <- rmd_summary("Continuous Delay-Differential")

  # Life History
  LH_section <- c(rmd_LAA(age = "1:info$LH$maxage", header = "## Life History\n"), 
                  rmd_WAA(age = "1:info$LH$maxage"), rmd_LW(),
                  rmd_mat(age = "1:info$LH$maxage", mat = "ifelse(1:info$LH$maxage < info$data$k, 0, 1)", 
                          fig.cap = "Assumed knife-edge maturity at age corresponding to length of 50% maturity."))
  
  # Data section
  if (any(Assessment@obj$env$data$MW_hist > 0, na.rm = TRUE)) {
    data_MW <- rmd_data_MW()
  } else {
    data_MW <- ""
  }
  data_section <- c(rmd_data_timeseries("Catch", header = "## Data\n"),
                    rmd_data_timeseries("Index", is_matrix = is.matrix(Assessment@Obs_Index), nsets = ncol(Assessment@Obs_Index)),
                    data_MW)

  # Assessment
  #### Pars and Fit
  assess_fit <- c(rmd_R0(header = "## Assessment {.tabset}\n### Estimates and Model Fit\n"), rmd_h(), rmd_M_prior(),
                  rmd_sel(age = "1:info$LH$maxage", sel = "ifelse(1:info$LH$maxage < info$data$k, 0, 1)", 
                          fig.cap = "Knife-edge selectivity set to the age corresponding to the length of 50% maturity."),
                  rmd_assess_fit_series(nsets = ncol(Assessment@Index)), rmd_assess_fit("Catch", "catch", match = TRUE))
  
  if (any(Assessment@obj$env$data$MW_hist > 0, na.rm = TRUE)) {
    fit_MW <- rmd_assess_fit_MW()
  } else {
    fit_MW <- ""
  }
  assess_fit <- c(assess_fit, fit_MW)
  
  if (state_space) {
    assess_fit2 <- c(rmd_residual("Dev", fig.cap = "Time series of recruitment deviations.", label = Assessment@Dev_type),
                     rmd_residual("Dev", "SE_Dev", fig.cap = "Time series of recruitment deviations with 95% confidence intervals.",
                                  label = Assessment@Dev_type, conv_check = TRUE))
    assess_fit <- c(assess_fit, assess_fit2)
  }

  #### Time Series
  ts_output <- c(rmd_F(header = "### Time Series Output\n"), rmd_F_FMSY(), rmd_SSB(), rmd_SSB_SSBMSY(),
                 rmd_SSB_SSB0(), rmd_dynamic_SSB0("TMB_report$dynamic_SSB0"),
                 rmd_Kobe("SSB_SSBMSY", xlab = "expression(SSB/SSB[MSY])"), rmd_R(),
                 rmd_N())

  #### Productivity
  SR_calc <- c("SSB_SR <- SSB[1:info$data$ny]",
               "R_SR <- R_pred(SSB_SR, TMB_report$h, TMB_report$R0, TMB_report$B0, info$data$SR_type)",
               "Rest <- R[1:info$data$ny + info$data$k]")
  productivity <- c(rmd_SR(header = "### Productivity\n\n\n", SR_calc = SR_calc),
                    rmd_SR(fig.cap = "Stock-recruit relationship (trajectory plot).", trajectory = TRUE),
                    rmd_yield_F("cDD"), rmd_yield_depletion("cDD"), rmd_sp(), rmd_SPR(), rmd_YPR())

  return(c(ss, LH_section, data_section, assess_fit, ts_output, productivity))
}


profile_likelihood_cDD <- function(Assessment, ...) {
  dots <- list(...)
  if (!"R0" %in% names(dots) && !"h" %in% names(dots)) stop("Sequence of neither R0 nor h was not found. See help file.")
  if (!is.null(dots$R0)) R0 <- dots$R0 else {
    R0 <- Assessment@R0
    profile_par <- "h"
  }
  if (!is.null(dots$h)) h <- dots$h else {
    h <- Assessment@h
    profile_par <- "R0"
  }

  map <- Assessment@obj$env$map
  params <- Assessment@info$params

  profile_grid <- expand.grid(R0 = R0, h = h)
  joint_profile <- !exists("profile_par", inherits = FALSE)

  profile_fn <- function(i, Assessment, params, map) {
    params$R0x <- log(profile_grid[i, 1] * Assessment@obj$env$data$rescale)

    if (Assessment@info$data$SR_type == "BH") {
      params$transformed_h <- logit((profile_grid[i, 2] - 0.2)/0.8)
    } else {
      params$transformed_h <- log(profile_grid[i, 2] - 0.2)
    }

    if (length(Assessment@opt$par) == 1) { # R0 is the only estimated parameter
      if (!joint_profile && profile_par == "R0") {
        nll <- Assessment@obj$fn(params$R0x)
      } else {

        obj2 <- MakeADFun(data = Assessment@info$data, parameters = params, map = map, random = Assessment@obj$env$random,
                          DLL = "SAMtool", silent = TRUE)

        if (joint_profile) {
          nll <- obj2$fn(params$R0x)
        } else { # Profile h
          opt2 <- optimize_TMB_model(obj2, Assessment@info$control, do_sd = FALSE)[[1]]
          if (!is.character(opt2)) nll <- opt2$objective else nll <- NA
        }
      }
    } else if (length(Assessment@opt$par) == 2) {

      if (all(names(Assessment@opt$par) == c("R0x", "transformed_h"))) {
        if (joint_profile) {
          nll <- Assessment@obj$fn(c(params$R0x, params$transformed_h))
        } else {
          if (profile_par == "R0") map$R0x <- factor(NA) else map$transformed_h <- factor(NA)
          obj2 <- MakeADFun(data = Assessment@info$data, parameters = params, map = map, random = Assessment@obj$env$random,
                            DLL = "SAMtool", silent = TRUE)
          opt2 <- optimize_TMB_model(obj2, Assessment@info$control, do_sd = FALSE)[[1]]
          if (!is.character(opt2)) nll <- opt2$objective else nll <- NA
        }
      } else { # R0, F
        if (joint_profile || profile_par == "R0") map$R0x <- factor(NA)
        obj2 <- MakeADFun(data = Assessment@info$data, parameters = params, map = map, random = Assessment@obj$env$random,
                          DLL = "SAMtool", silent = TRUE)
        opt2 <- optimize_TMB_model(obj2, Assessment@info$control, do_sd = FALSE)[[1]]
        if (!is.character(opt2)) nll <- opt2$objective else nll <- NA
      }

    } else { # more than 2 parameters

      if (joint_profile) {
        map$R0x <- factor(NA)
        if (!"transformed_h" %in% names(Assessment@opt$par)) map$transformed_h <- factor(NA)
      } else {
        if (profile_par == "R0") map$R0x <- factor(NA) else map$transformed_h <- factor(NA)
      }

      obj2 <- MakeADFun(data = Assessment@info$data, parameters = params, map = map, random = Assessment@obj$env$random, DLL = "SAMtool", silent = TRUE)
      opt2 <- optimize_TMB_model(obj2, Assessment@info$control, do_sd = FALSE)[[1]]
      if (!is.character(opt2)) nll <- opt2$objective else nll <- NA_real_
    }
    if (!exists("nll", inherits = FALSE)) nll <- NA_real_
    return(nll)
  }
  nll <- pblapply(1:nrow(profile_grid), profile_fn, 
                  Assessment = Assessment, params = params, map = map,
                  cl = if (snowfall::sfIsRunning()) snowfall::sfGetCluster() else NULL)
  profile_grid$nll <- do.call(c, nll) - Assessment@opt$objective

  if (joint_profile) {
    pars <- c("R0", "h")
    MLE <- vapply(pars, function(x, y) slot(y, x), y = Assessment, numeric(1))
  } else {
    pars <- profile_par
    MLE <- slot(Assessment, pars)
  }

  output <- new("prof", Model = Assessment@Model, Name = Assessment@Name, Par = pars, MLE = MLE, grid = profile_grid)
  return(output)
}


retrospective_cDD <- function(Assessment, nyr, state_space = FALSE) {
  assign_Assessment_slots(Assessment)
  ny <- info$data$ny
  k <- info$data$k

  Year <- info$Year
  moreRecruitYears <- max(Year) + 1:k
  Year <- c(Year, moreRecruitYears)

  # Array dimension: Retroyr, Year, ts
  # ts includes: F, F/FMSY, B, B/BMSY, B/B0, R, VB
  retro_ts <- array(NA, dim = c(nyr+1, ny+k, 7))
  TS_var <- c("F", "F_FMSY", "B", "B_BMSY", "B_B0", "R", "VB")
  dimnames(retro_ts) <- list(Peel = 0:nyr, Year = Year, Var = TS_var)

  retro_est <- array(NA, dim = c(nyr+1, length(SD$par.fixed[names(SD$par.fixed) != "log_rec_dev"]), 2))
  dimnames(retro_est) <- list(Peel = 0:nyr, Var = names(SD$par.fixed)[names(SD$par.fixed) != "log_rec_dev"],
                              Value = c("Estimate", "Std. Error"))

  lapply_fn <- function(i, info, obj, state_space) {
    ny_ret <- info$data$ny - i
    info$data$ny <- ny_ret
    info$data$C_hist <- info$data$C_hist[1:ny_ret]
    info$data$I_hist <- info$data$I_hist[1:ny_ret, , drop = FALSE]
    info$data$I_sd <- info$data$I_sd[1:ny_ret, , drop = FALSE]
    info$data$MW_hist <- info$data$MW_hist[1:ny_ret]

    if (state_space) info$params$log_rec_dev <- rep(0, ny_ret)

    obj2 <- MakeADFun(data = info$data, parameters = info$params, map = obj$env$map, random = obj$env$random,
                      inner.control = info$inner.control, DLL = "SAMtool", silent = TRUE)
    mod <- optimize_TMB_model(obj2, info$control)
    opt2 <- mod[[1]]
    SD <- mod[[2]]

    if (!is.character(opt2)) {
      report <- obj2$report(obj2$env$last.par.best)
      ref_pt <- ref_pt_cDD(info$data, report$Arec, report$Brec, report$M)
      report <- c(report, ref_pt)

      FMort <- c(report$F, rep(NA, k + i))
      F_FMSY <- FMort/report$FMSY
      B <- c(report$B, rep(NA, k - 1 + i))
      B_BMSY <- B/report$BMSY
      B_B0 <- B/report$B0
      R <- c(report$R, rep(NA, i))
      VB <- B

      retro_ts[i+1, , ] <<- cbind(FMort, F_FMSY, B, B_BMSY, B_B0, R, VB)

      sumry <- summary(SD, "fixed")
      sumry <- sumry[rownames(sumry) != "log_rec_dev", drop = FALSE]
      retro_est[i+1, , ] <<- sumry

      return(SD$pdHess)
    }
    return(FALSE)
  }

  conv <- vapply(0:nyr, lapply_fn, logical(1), info = info, obj = obj, state_space = state_space)
  if (any(!conv)) warning("Peels that did not converge: ", paste0(which(!conv) - 1, collapse = " "))

  retro <- new("retro", Model = Assessment@Model, Name = Assessment@Name, TS_var = TS_var, TS = retro_ts,
               Est_var = dimnames(retro_est)[[2]], Est = retro_est)
  attr(retro, "TS_lab") <- c("Fishing mortality", expression(F/F[MSY]), "Biomass", expression(B/B[MSY]), expression(B/B[0]),
                             "Recruitment", "Vulnerable biomass")

  return(retro)
}



summary_cDD_SS <- function(Assessment) summary_cDD(Assessment, TRUE)

rmd_cDD_SS <- function(Assessment, ...) rmd_cDD(Assessment, TRUE, ...)

profile_likelihood_cDD_SS <- profile_likelihood_cDD

retrospective_cDD_SS <- function(Assessment, nyr) retrospective_cDD(Assessment, nyr, TRUE)


plot_yield_cDD <- function(data, report, fmsy, msy, xaxis = c("F", "Biomass", "Depletion")) {
  xaxis <- match.arg(xaxis)
  if (xaxis == "F") F.vector <- seq(0, 2.5 * fmsy, length.out = 1e2) else F.vector <- seq(0, 5 * fmsy, length.out = 1e2)
  
  yield <- lapply(F.vector, yield_fn_cDD, M = report$M, Kappa = data$Kappa, 
                  Winf = data$Winf, wk = data$wk, SR = data$SR_type, 
                  Arec = report$Arec, Brec = report$Brec, opt = FALSE)
  Biomass <- vapply(yield, getElement, numeric(1), "B")
  Yield <- vapply(yield, getElement, numeric(1), "Yield")
  R <- vapply(yield, getElement, numeric(1), "R")
  ind <- R >= 0
  
  if (missing(fmsy)) fmsy <- F.vector[which.max(Yield)[1]]
  if (missing(msy)) msy <- max(Yield)

  if (xaxis == "F") {
    plot(F.vector[ind], Yield[ind], typ = 'l', xlab = "Fishing Mortality", ylab = "Equilibrium yield")
    segments(x0 = fmsy, y0 = 0, y1 = msy, lty = 2)
    segments(x0 = 0, y0 = msy, x1 = fmsy, lty = 2)
    abline(h = 0, col = 'grey')
  }

  if (xaxis == "Biomass") {
    plot(Biomass[ind], Yield[ind], typ = 'l', xlab = "Biomass", ylab = "Equilibrium yield")
    segments(x0 = report$BMSY, y0 = 0, y1 = msy, lty = 2)
    segments(x0 = 0, y0 = msy, x1 = report$BMSY, lty = 2)
    abline(h = 0, col = 'grey')
  }

  if (xaxis == "Depletion") {
    plot(Biomass[ind]/report$B0, Yield[ind], typ = 'l', xlab = expression(B/B[0]), ylab = "Equilibrium yield")
    segments(x0 = report$BMSY/report$B0, y0 = 0, y1 = msy, lty = 2)
    segments(x0 = 0, y0 = msy, x1 = report$BMSY/report$B0, lty = 2)
    abline(h = 0, col = 'grey')
  }
  invisible(data.frame(F = F.vector[ind], Yield = Yield[ind], B = Biomass[ind], B_B0 = Biomass[ind]/report$B0))
}
