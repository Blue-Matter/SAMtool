
summary_DD_TMB <- function(Assessment, state_space = FALSE) {
  assign_Assessment_slots()

  if(conv) current_status <- c(U_UMSY[length(U_UMSY)], B_BMSY[length(B_BMSY)], B_B0[length(B_B0)])
  else current_status <- c(NA, NA, B_B0[length(B_B0)])
  current_status <- data.frame(Value = current_status)
  rownames(current_status) <- c("U/UMSY", "B/BMSY", "B/B0")

  Value <- c(unlist(info$data[c(2,3,4,6,7)]))
  Description <- c("Unfished survival = exp(-M)", "alpha = Winf * (1-rho)",
                  "rho = (W_k+2 - Winf)/(W_k+1 - Winf)",
                  "Age of knife-edge selectivity",
                  "Weight at age k")
  rownam <- c("S0", "alpha", "rho", "k", "w_k")
  if(Assessment@obj$env$data$condition == "effort" && "log_omega" %in% names(obj$env$map)) {
    Value <- c(Value, TMB_report$omega)
    Description <- c(Description, "Catch SD (log-space)")
    rownam <- c(rownam, "omega")
  }
  if(state_space && "log_tau" %in% names(obj$env$map)) {
    Value <- c(Value, TMB_report$tau)
    Description <- c(Description, "log-Recruitment deviation SD")
    rownam <- c(rownam, "tau")
  }
  if("transformed_h" %in% names(obj$env$map)) {
    Value <- c(Value, h)
    Description <- c(Description, "Stock-recruit steepness")
    rownam <- c(rownam, "h")
  }
  input_parameters <- data.frame(Value = Value, Description = Description, stringsAsFactors = FALSE)
  rownames(input_parameters) <- rownam

  if(conv) derived <- c(B0, N0, MSY, UMSY, BMSY, BMSY/B0)
  else derived <- rep(NA, 6)
  derived <- data.frame(Value = derived,
                        Description = c("Unfished biomass", "Unfished abundance", "Maximum sustainable yield (MSY)",
                                        "Harvest rate at MSY", "Biomass at MSY", "Depletion at MSY"),
                        stringsAsFactors = FALSE)
  rownames(derived) <- c("B0", "N0", "MSY", "UMSY", "BMSY", "BMSY/B0")

  model_estimates <- sdreport_int(SD)
  if(!is.character(model_estimates)) {
    rownames(model_estimates)[rownames(model_estimates) == "log_rec_dev"] <- paste0("log_rec_dev_", names(Dev))
  }

  model_name <- "Delay-Difference"
  if(state_space) model_name <- paste(model_name, "(State-Space)")
  output <- list(model = model_name,
                 current_status = current_status, input_parameters = input_parameters,
                 derived_quantities = derived, model_estimates = model_estimates,
                 log_likelihood = matrix(NLL, ncol = 1, dimnames = list(names(NLL), "Neg.LL")))
  return(output)
}


rmd_DD_TMB <- function(Assessment, state_space = FALSE, ...) {
  if(state_space) {
    ss <- rmd_summary("Delay-Difference (State-Space)")
  } else ss <- rmd_summary("Delay-Difference")

  # Life History
  age <- 1:Assessment@info$LH$maxage
  k <- Assessment@info$data$k
  mat <- ifelse(age < k, 0, 1)
  LH_section <- c(rmd_LAA(age, Assessment@info$LH$LAA, header = "## Life History\n"), rmd_WAA(age, Assessment@info$LH$WAA),
                  rmd_LW(Assessment@info$LH$LAA, Assessment@info$LH$WAA),
                  rmd_mat(age, mat, fig.cap = "Assumed knife-edge maturity at age corresponding to length of 50% maturity."))

  # Data section
  data_section <- c(rmd_data_timeseries("Catch", header = "## Data\n"),
                    rmd_data_timeseries("Index", is_matrix = is.matrix(Assessment@Obs_Index), nsets = ncol(Assessment@Obs_Index)))

  # Assessment
  #### Pars and Fit
  assess_all <- c(rmd_R0(header = "## Assessment {.tabset}\n### Estimates and Model Fit\n"), rmd_h(),
                  rmd_sel(age, mat, fig.cap = "Knife-edge selectivity set to the age corresponding to the length of 50% maturity."))

  if(Assessment@obj$env$data$condition == "effort") {
    assess_data <- c(rmd_assess_fit("Catch", "catch"), rmd_assess_resid("Catch"), rmd_assess_qq("Catch", "catch"))
  } else {
    assess_data <- rmd_assess_fit_series(nsets = ncol(Assessment@Index))
  }
  assess_fit <- c(assess_all, assess_data)

  if(state_space) {
    assess_fit2 <- c(rmd_residual("Dev", fig.cap = "Time series of recruitment deviations.", label = Assessment@Dev_type),
                     rmd_residual("Dev", "SE_Dev", fig.cap = "Time series of recruitment deviations with 95% confidence intervals.",
                                  label = Assessment@Dev_type, conv_check = TRUE))
    assess_fit <- c(assess_fit, assess_fit2)
  }

  #### Time Series
  ts_output <- c(rmd_U(header = "### Time Series Output\n"), rmd_U_UMSY(), rmd_SSB(), rmd_SSB_SSBMSY(),
                 rmd_SSB_SSB0(), rmd_Kobe("SSB_SSBMSY", "U_UMSY", xlab = "expression(SSB/SSB[MSY])", ylab = "expression(U/U[MSY])"),
                 rmd_R(), rmd_N())

  #### Productivity
  ny <- Assessment@info$data$ny
  SSB <- Assessment@SSB[1:ny]
  Arec <- Assessment@TMB_report$Arec
  Brec <- Assessment@TMB_report$Brec
  if(Assessment@info$data$SR_type == "BH") expectedR <- Arec * SSB / (1 + Brec * SSB) else {
    expectedR <- Arec * SSB * exp(-Brec * SSB)
  }

  first_recruit_year <- k + 1
  last_recruit_year <- length(Assessment@info$Year) + k
  ind_recruit <- first_recruit_year:last_recruit_year
  rec_dev <- Assessment@R[ind_recruit]

  productivity <- c(rmd_SR(SSB, expectedR, rec_dev, header = "### Productivity\n\n\n"),
                    rmd_SR(SSB, expectedR, rec_dev, fig.cap = "Stock-recruit relationship (trajectory plot).", trajectory = TRUE),
                    rmd_yield_U("DD"), rmd_yield_depletion("DD"), rmd_sp())

  return(c(ss, LH_section, data_section, assess_fit, ts_output, productivity))

}



profile_likelihood_DD_TMB <- function(Assessment, ...) {
  dots <- list(...)
  if(!"R0" %in% names(dots) && !"h" %in% names(dots)) stop("Sequence of neither R0 nor h was not found. See help file.")
  if(!is.null(dots$R0)) R0 <- dots$R0 else {
    R0 <- Assessment@R0
    profile_par <- "h"
  }
  if(!is.null(dots$h)) h <- dots$h else {
    h <- Assessment@h
    profile_par <- "R0"
  }

  map <- Assessment@obj$env$map
  params <- Assessment@info$params

  profile_grid <- expand.grid(R0 = R0, h = h)
  joint_profile <- !exists("profile_par")

  profile_fn <- function(i, Assessment, params, map) {
    params$R0x <- log(profile_grid[i, 1] * Assessment@obj$env$data$rescale)
    if(Assessment@info$data$SR_type == "BH") {
      params$transformed_h <- logit((profile_grid[i, 2] - 0.2)/0.8)
    } else {
      params$transformed_h <- log(profile_grid[i, 2] - 0.2)
    }

    if(joint_profile) {
      map$R0x <- map$transformed_h <- factor(NA)
    } else {
      if(profile_par == "R0") map$R0x <- factor(NA) else map$transformed_h <- factor(NA)
    }
    obj2 <- MakeADFun(data = Assessment@info$data, parameters = params, map = map, random = Assessment@obj$env$random,
                      DLL = "MSEtool", silent = TRUE)
    opt2 <- optimize_TMB_model(obj2, Assessment@info$control)[[1]]
    if(!is.character(opt2)) nll <- opt2$objective else nll <- NA
    return(nll)
  }
  nll <- vapply(1:nrow(profile_grid), profile_fn, numeric(1), Assessment = Assessment, params = params, map = map) - Assessment@opt$objective
  profile_grid$nll <- nll

  if(joint_profile) {
    pars <- c("R0", "h")
    MLE <- vapply(pars, function(x, y) slot(y, x), y = Assessment, numeric(1))
  } else {
    pars <- profile_par
    MLE <- slot(Assessment, pars)
  }

  output <- new("prof", Model = Assessment@Model, Name = Assessment@Name, Par = pars, MLE = MLE, grid = profile_grid)
  return(output)
}


retrospective_DD_TMB <- function(Assessment, nyr, state_space = FALSE) {
  assign_Assessment_slots(Assessment)
  ny <- info$data$ny
  k <- info$data$k

  Year <- info$Year
  moreRecruitYears <- max(Year) + 1:k
  Year <- c(Year, moreRecruitYears)

  # Array dimension: Retroyr, Year, ts
  # ts includes: U, U/UMSY, B, B/BMSY, B/B0, R, VB
  retro_ts <- array(NA, dim = c(nyr+1, ny+k, 7))
  TS_var <- c("U", "U_UMSY", "B", "B_BMSY", "B_B0", "R", "VB")
  dimnames(retro_ts) <- list(Peel = 0:nyr, Year = Year, Var = TS_var)

  retro_est <- array(NA, dim = c(nyr+1, length(SD$par.fixed[names(SD$par.fixed) != "log_rec_dev"]), 2))
  dimnames(retro_est) <- list(Peel = 0:nyr, Var = names(SD$par.fixed)[names(SD$par.fixed) != "log_rec_dev"],
                              Value = c("Estimate", "Std. Error"))

  lapply_fn <- function(i, info, obj, state_space) {
    ny_ret <- ny - i
    info$data$ny <- ny_ret
    info$data$C_hist <- info$data$C_hist[1:ny_ret]
    info$data$E_hist <- info$data$E_hist[1:ny_ret]
    info$data$I_hist <- info$data$I_hist[1:ny_ret, , drop = FALSE]
    info$data$I_sd <- info$data$I_sd[1:ny_ret, , drop = FALSE]

    if(state_space) info$params$log_rec_dev <- rep(0, ny_ret - k)

    obj2 <- MakeADFun(data = info$data, parameters = info$params, random = obj$env$random, map = obj$env$map,
                      inner.control = info$inner.control, DLL = "MSEtool", silent = TRUE)
    mod <- optimize_TMB_model(obj2, info$control)
    opt2 <- mod[[1]]
    SD <- mod[[2]]

    if(!is.character(opt2) && !is.character(SD)) {
      report <- obj2$report(obj2$env$last.par.best)
      ref_pt <- get_MSY_DD(info$data, report$Arec, report$Brec)
      report <- c(report, ref_pt)

      U <- c(report$U, rep(NA, k + i))
      U_UMSY <- U/report$UMSY
      B <- c(report$B, rep(NA, k - 1 + i))
      B_BMSY <- B/report$BMSY
      B_B0 <- B/B0
      R <- c(report$R, rep(NA, i))
      VB <- B

      retro_ts[i+1, , ] <<- cbind(U, U_UMSY, B, B_BMSY, B_B0, R, VB)

      sumry <- summary(SD, "fixed")
      sumry <- sumry[rownames(sumry) != "log_rec_dev", drop = FALSE]
      retro_est[i+1, , ] <<- sumry

      return(SD$pdHess)
    }
    return(FALSE)
  }

  conv <- vapply(0:nyr, lapply_fn, logical(1), info = info, obj = obj, state_space = state_space)
  if(any(!conv)) warning("Peels that did not converge: ", paste0(which(!conv) - 1, collapse = " "))

  retro <- new("retro", Model = Assessment@Model, Name = Assessment@Name, TS_var = TS_var, TS = retro_ts,
               Est_var = dimnames(retro_est)[[2]], Est = retro_est)
  attr(retro, "TS_lab") <- c("Harvest rate", expression(U/U[MSY]), "Biomass", expression(B/B[MSY]), expression(B/B[0]),
                             "Recruitment", "Vulnerable biomass")

  return(retro)
}



summary_DD_SS <- function(Assessment) summary_DD_TMB(Assessment, TRUE)

rmd_DD_SS <- function(Assessment, ...) rmd_DD_TMB(Assessment, TRUE, ...)

profile_likelihood_DD_SS <- profile_likelihood_DD_TMB

retrospective_DD_SS <- function(Assessment, nyr) retrospective_DD_TMB(Assessment, nyr, TRUE)


plot_yield_DD <- function(data, report, umsy, msy, xaxis = c("U", "Biomass", "Depletion")) {
  xaxis <- match.arg(xaxis)
  u.vector <- seq(0, 1, 0.01)
  S0 <- data$S0
  Alpha <- data$Alpha
  wk <- data$wk
  Rho <- data$Rho
  SR_type  <- data$SR_type

  Arec <- report$Arec
  Brec <- report$Brec
  BMSY <- report$BMSY

  Surv <- S0 * (1 - u.vector)

  BPR <- (Surv * Alpha/(1 - Surv) + wk)/(1 - Rho * Surv)
  if(SR_type == "BH") R <- (Arec * BPR - 1)/(Brec * BPR)
  if(SR_type == "Ricker") R <- log(Arec * BPR)/(Brec * BPR)

  Biomass <- BPR * R
  Yield <- u.vector * BPR * R
  ind <- R >= 0

  if(xaxis == "U") {
    plot(u.vector[ind], Yield[ind], typ = 'l', xlab = "Exploitation rate (U)",
         ylab = "Equilibrium yield")
    segments(x0 = umsy, y0 = 0, y1 = msy, lty = 2)
    segments(x0 = 0, y0 = msy, x1 = umsy, lty = 2)
    abline(h = 0, col = 'grey')
  }

  if(xaxis == "Biomass") {
    plot(Biomass[ind], Yield[ind], typ = 'l', xlab = "Biomass",
         ylab = "Equilibrium yield")
    segments(x0 = BMSY, y0 = 0, y1 = msy, lty = 2)
    segments(x0 = 0, y0 = msy, x1 = BMSY, lty = 2)
    abline(h = 0, col = 'grey')
  }

  if(xaxis == "Depletion") {
    plot(Biomass[ind]/report$B0, Yield[ind], typ = 'l',
         xlab = expression(B/B[0]), ylab = "Equilibrium yield")
    segments(x0 = BMSY/report$B0, y0 = 0, y1 = msy, lty = 2)
    segments(x0 = 0, y0 = msy, x1 = BMSY/report$B0, lty = 2)
    abline(h = 0, col = 'grey')
  }
  invisible(data.frame(U = u.vector[ind], Yield = Yield[ind], B = Biomass[ind], B_B0 = Biomass[ind]/report$B0))
}


