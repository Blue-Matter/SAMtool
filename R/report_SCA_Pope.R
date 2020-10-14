
summary_SCA_Pope <- function(Assessment) {
  assign_Assessment_slots(Assessment)

  if(conv) current_status <- c(U_UMSY[length(U_UMSY)], SSB_SSBMSY[length(SSB_SSBMSY)], SSB_SSB0[length(SSB_SSB0)])
  else current_status <- c(NA, NA, SSB_SSB0[length(SSB_SSB0)])
  current_status <- data.frame(Value = current_status)
  rownames(current_status) <- c("U/UMSY", "SSB/SSBMSY", "SSB/SSB0")

  Value <- c(h, info$data$M[1], info$data$max_age, info$LH$Linf, info$LH$K, info$LH$t0,
             info$LH$a * info$LH$Linf ^ info$LH$b, info$LH$A50, info$LH$A95)
  Description = c("Stock-recruit steepness", "Natural mortality", "Maximum age (plus-group)", "Asymptotic length", "Growth coefficient",
                  "Age at length-zero", "Asymptotic weight", "Age of 50% maturity", "Age of 95% maturity")
  rownam <- c("h", "M", "maxage", "Linf", "K", "t0", "Winf", "A50", "A95")
  input_parameters <- data.frame(Value = Value, Description = Description, stringsAsFactors = FALSE)
  rownames(input_parameters) <- rownam
  if(!"transformed_h" %in% names(obj$env$map)) input_parameters <- input_parameters[-1, ]

  if(conv) Value <- c(VB0, SSB0, MSY, UMSY, VBMSY, SSBMSY, SSBMSY/SSB0)
  else Value <- rep(NA, 7)

  Description <- c("Unfished vulnerable biomass",
                  "Unfished spawning stock biomass (SSB)", "Maximum sustainable yield (MSY)", "Harvest rate at MSY",
                  "Vulnerable biomass at MSY", "SSB at MSY", "Spawning depletion at MSY")
  derived <- data.frame(Value = Value, Description = Description, stringsAsFactors = FALSE)
  rownames(derived) <- c("VB0", "SSB0", "MSY", "UMSY", "VBMSY", "SSBMSY", "SSBMSY/SSB0")

  if(!is.character(model_estimates)) {
    rownames(model_estimates)[rownames(model_estimates) == "log_rec_dev"] <- paste0("log_rec_dev_", names(FMort)[as.logical(obj$env$data$est_rec_dev)])
  }

  output <- list(model = "Statistical Catch-at-Age (SCA_Pope)",
                 current_status = current_status, input_parameters = input_parameters,
                 derived_quantities = derived, model_estimates = model_estimates,
                 log_likelihood = matrix(NLL, ncol = 1, dimnames = list(names(NLL), "Neg.LL")))
  return(output)
}

rmd_SCA_Pope <- function(Assessment, ...) {
  ss <- rmd_summary("Statistical Catch-at-Age (SCA_Pope)")

  # Life History
  age <- 1:Assessment@info$data$max_age
  LH_section <- c(rmd_LAA(age, Assessment@info$LH$LAA, header = "## Life History\n"), rmd_WAA(age, Assessment@info$LH$WAA),
                  rmd_LW(Assessment@info$LH$LAA, Assessment@info$LH$WAA),
                  rmd_mat(age, Assessment@info$data$mat,
                          fig.cap = "Maturity at age. Length-based maturity parameters were converted to the corresponding ages."))
  # Data section
  data_section <- c(rmd_data_timeseries("Catch", header = "## Data\n"), rmd_data_timeseries("Index"),
                    rmd_data_age_comps("bubble"), rmd_data_age_comps("annual"))

  # Assessment
  #### Pars and Fit
  assess_fit <- c(rmd_R0(header = "## Assessment {.tabset}\n### Estimates and Model Fit\n"), rmd_h(),
                  rmd_sel(age, Assessment@Selectivity[nrow(Assessment@Selectivity), ], fig.cap = "Estimated selectivity at age."),
                  rmd_assess_fit("Index", "index"), rmd_assess_resid("Index"), rmd_assess_qq("Index", "index"),
                  rmd_fit_age_comps("bubble"), rmd_fit_age_comps("annual"),
                  rmd_residual("Dev", fig.cap = "Time series of recruitment deviations.", label = Assessment@Dev_type,
                               blue = any(as.numeric(names(Assessment@Dev)) < Assessment@info$Year[1])),
                  rmd_residual("Dev", "SE_Dev", fig.cap = "Time series of recruitment deviations with 95% confidence intervals.",
                               label = Assessment@Dev_type, conv_check = TRUE, blue = any(as.numeric(names(Assessment@Dev)) < Assessment@info$Year[1])))

  #### Time Series
  ts_output <- c(rmd_U(header = "### Time Series Output\n"), rmd_U_UMSY(), rmd_SSB(), rmd_SSB_SSBMSY(),
                 rmd_SSB_SSB0(), rmd_Kobe("SSB_SSBMSY", "U_UMSY", xlab = "expression(SSB/SSB[MSY])", ylab = "expression(U/U[MSY])"), rmd_R(),
                 rmd_N(), rmd_N_at_age(), rmd_C_at_age(), rmd_C_mean_age())

  # Productivity
  Arec <- Assessment@TMB_report$Arec
  Brec <- Assessment@TMB_report$Brec
  SSB <- Assessment@SSB[1:(length(Assessment@SSB)-1)]

  SR <- Assessment@info$data$SR_type
  if(SR == "BH") expectedR <- Arec * SSB / (1 + Brec * SSB) else {
    expectedR <- Arec * SSB * exp(-Brec * SSB)
  }
  estR <- Assessment@R[as.numeric(names(Assessment@R)) > Assessment@info$Year[1]]

  productivity <- c(rmd_SR(SSB, expectedR, estR, header = "### Productivity\n\n\n"),
                    rmd_SR(SSB, expectedR, estR, fig.cap = "Stock-recruit relationship (trajectory plot).", trajectory = TRUE),
                    rmd_yield_U("SCA_Pope"), rmd_yield_depletion("SCA_Pope"), rmd_sp())

  return(c(ss, LH_section, data_section, assess_fit, ts_output, productivity))
}


profile_likelihood_SCA_Pope <- profile_likelihood_SCA


retrospective_SCA_Pope <- function(Assessment, nyr) {
  assign_Assessment_slots(Assessment)
  n_y <- info$data$n_y

  Year <- c(info$Year, max(info$Year) + 1)

  # Array dimension: Retroyr, Year, ts
  # ts includes: Calendar U, U_MSY, B, B/BMSY, B/B0, R, VB
  retro_ts <- array(NA, dim = c(nyr+1, n_y + 1, 7))
  TS_var <- c("U", "U_UMSY", "SSB", "SSB_SSBMSY", "SSB_SSB0", "R", "VB")
  dimnames(retro_ts) <- list(Peel = 0:nyr, Year = Year, Var = TS_var)

  SD_nondev <- summary(SD)[rownames(summary(SD)) != "log_rec_dev" & rownames(summary(SD)) != "log_early_rec_dev", ]
  retro_est <- array(NA, dim = c(nyr+1, dim(SD_nondev)))
  dimnames(retro_est) <- list(Peel = 0:nyr, Var = rownames(SD_nondev), Value = c("Estimate", "Std. Error"))

  lapply_fn <- function(i, info, obj) {
    n_y_ret <- n_y - i
    info$data$n_y <- n_y_ret
    info$data$C_hist <- info$data$C_hist[1:n_y_ret]
    info$data$I_hist <- info$data$I_hist[1:n_y_ret]
    info$data$CAA_hist <- info$data$CAA_hist[1:n_y_ret, ]
    info$data$CAA_n <- info$data$CAA_n[1:n_y_ret]
    info$data$est_rec_dev <- info$data$est_rec_dev[1:n_y_ret]

    info$params$log_rec_dev <- rep(0, n_y_ret)

    map <- obj$env$map
    if(any(names(map) == "log_rec_dev")) {
      new_map <- as.numeric(map$log_rec_dev) - i
      map$log_rec_dev <- factor(new_map[new_map > 0])
    }

    obj2 <- MakeADFun(data = info$data, parameters = info$params, map = map, random = obj$env$random,
                      inner.control = info$inner.control, DLL = "MSEtool", silent = TRUE)
    mod <- optimize_TMB_model(obj2, info$control)
    opt2 <- mod[[1]]
    SD <- mod[[2]]

    if(!is.character(opt2) && !is.character(SD)) {
      report <- obj2$report(obj2$env$last.par.best)
      ref_pt <- SCA_Pope_MSY_calc(Arec = report$Arec, Brec = report$Brec, M = info$data$M, weight = info$data$weight, mat = info$data$mat,
                                  vul = report$vul, SR = info$data$SR_type)

      report <- c(report, ref_pt)

      U <- c(report$U, rep(NA, i + 1))
      U_UMSY <- U/report$UMSY
      SSB <- c(report$E, rep(NA, i))
      SSB_SSBMSY <- SSB/report$EMSY
      SSB_SSB0 <- SSB/report$E0
      R <- c(report$R, rep(NA, i))
      VB <- c(report$E, rep(NA, i))
      #log_rec_dev <- c(report$log_rec_dev, rep(NA, i + 1))

      retro_ts[i+1, , ] <<- cbind(U, U_UMSY, SSB, SSB_SSBMSY, SSB_SSB0, R, VB)
      retro_est[i+1, , ] <<- summary(SD)[rownames(summary(SD)) != "log_rec_dev" & rownames(summary(SD)) != "log_early_rec_dev", ]

      return(SD$pdHess)
    }
    return(FALSE)
  }

  conv <- vapply(0:nyr, lapply_fn, logical(1), info = info, obj = obj)
  if(any(!conv)) warning("Peels that did not converge: ", paste0(which(!conv) - 1, collapse = " "))

  retro <- new("retro", Model = Assessment@Model, Name = Assessment@Name, TS_var = TS_var, TS = retro_ts,
               Est_var = dimnames(retro_est)[[2]], Est = retro_est)
  attr(retro, "TS_lab") <- c("Harvest rate", expression(U/U[MSY]), "Spawning biomass", expression(SSB/SSB[MSY]), "Spawning depletion",
                             "Recruitment", "Vulnerable biomass")

  return(retro)
}


plot_yield_SCA_Pope <- function(data, report, umsy, msy, xaxis = c("U", "Biomass", "Depletion")) {
  xaxis <- match.arg(xaxis)
  if(xaxis == "U") u.vector = seq(0, max(1, 2.5 * umsy), length.out = 100) else {
    u.vector = seq(0, 1, length.out = 100)
  }

  M <- data$M
  mat <- data$mat
  weight <- data$weight
  maxage <- data$max_age
  vul <- report$vul

  BMSY <- report$EMSY
  B0 <- report$E0

  SR <- data$SR_type

  Arec <- report$Arec
  Brec <- report$Brec

  EPR <- Req <- NA
  solveMSY <- function(logit_U) {
    U <- ilogit(logit_U)
    surv <- exp(-M) * (1 - vul * U)
    NPR <- c(1, cumprod(surv[1:(maxage-1)]))
    NPR[maxage] <- NPR[maxage]/(1 - surv[maxage])
    EPR <<- sum(NPR * mat * weight)
    if(SR == "BH") Req <<- (Arec * EPR - 1)/(Brec * EPR)
    if(SR == "Ricker") Req <<- log(Arec * EPR)/(Brec * EPR)
    CPR <- vul * U * NPR * exp(-0.5 * M)
    Yield <- Req * sum(CPR * weight)
    return(-1 * Yield)
  }

  Biomass <- Yield <- R <- rep(NA, length(u.vector))
  for(i in 1:length(u.vector)) {
    Yield[i] <- -1 * solveMSY(logit(u.vector[i]))
    R[i] <- Req
    Biomass[i] <- EPR * Req
  }

  ind <- R >= 0

  if(xaxis == "U") {
    plot(u.vector[ind], Yield[ind], typ = 'l', xlab = "Harvest rate (U)",
         ylab = "Equilibrium yield")
    segments(x0 = umsy, y0 = 0, y1 = msy, lty = 2)
    segments(x0 = 0, y0 = msy, x1 = umsy, lty = 2)
    abline(h = 0, col = 'grey')
  }

  if(xaxis == "Biomass") {
    plot(Biomass[ind], Yield[ind], typ = 'l', xlab = "Spawning Stock Biomass",
         ylab = "Equilibrium yield")
    segments(x0 = BMSY, y0 = 0, y1 = msy, lty = 2)
    segments(x0 = 0, y0 = msy, x1 = BMSY, lty = 2)
    abline(h = 0, col = 'grey')
  }

  if(xaxis == "Depletion") {
    plot(Biomass[ind]/B0, Yield[ind], typ = 'l',
         xlab = expression(SSB/SSB[0]), ylab = "Equilibrium yield")
    segments(x0 = BMSY/B0, y0 = 0, y1 = msy, lty = 2)
    segments(x0 = 0, y0 = msy, x1 = BMSY/report$B0, lty = 2)
    abline(h = 0, col = 'grey')
  }
  invisible(data.frame(U = u.vector[ind], Yield = Yield[ind], B = Biomass[ind], B_B0 = Biomass[ind]/report$B0))
}


