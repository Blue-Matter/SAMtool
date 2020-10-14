
summary_SSS <- function(Assessment) {
  output <- summary_SCA_Pope(Assessment)
  output$model <- "Simple Stock Synthesis (SSS)"
  return(output)
}

rmd_SSS <- function(Assessment, ...) {
  ss <- rmd_summary("Simple Stock Synthesis (SSS)")

  # Life History
  age <- 1:Assessment@info$data$max_age
  LH_section <- c(rmd_LAA(age, Assessment@info$LH$LAA, header = "## Life History\n"), rmd_WAA(age, Assessment@info$LH$WAA),
                  rmd_LW(Assessment@info$LH$LAA, Assessment@info$LH$WAA),
                  rmd_mat(age, Assessment@info$data$mat,
                          fig.cap = "Maturity at age. Length-based maturity parameters were converted to the corresponding ages."),
                  rmd_sel(age, Assessment@Selectivity[nrow(Assessment@Selectivity), ], fig.cap = "Selectivity at age."))

  # Data section
  data_section <- rmd_data_timeseries("Catch", header = "## Data\n")

  # Assessment
  #### Pars and Fit
  assess_fit <- c(rmd_R0(header = "## Assessment {.tabset}\n### Estimates and Model Fit\n"),
                  rmd_assess_fit("Index", "index"))

  #### Time Series
  ts_output <- c(rmd_U(header = "### Time Series Output\n"), rmd_U_UMSY(), rmd_SSB(), rmd_SSB_SSBMSY(),
                 rmd_SSB_SSB0(), rmd_Kobe("SSB_SSBMSY", "U_UMSY", xlab = "expression(SSB/SSB[MSY])", ylab = "expression(U/U[MSY])"), rmd_R(),
                 rmd_N(), rmd_N_at_age())

  # Productivity
  SSB <- Assessment@SSB[1:(length(Assessment@SSB)-1)]
  expectedR <- estR <- Assessment@R[as.numeric(names(Assessment@R)) > Assessment@info$Year[1]]

  productivity <- c(rmd_SR(SSB, expectedR, estR, header = "### Productivity\n\n\n"),
                    rmd_SR(SSB, expectedR, estR, fig.cap = "Stock-recruit relationship (trajectory plot).", trajectory = TRUE),
                    rmd_yield_U("SCA_Pope"), rmd_yield_depletion("SCA_Pope"), rmd_sp())

  return(c(ss, LH_section, data_section, assess_fit, ts_output, productivity))
}


profile_likelihood_SSS <- function(Assessment, ...) {
  dots <- list(...)
  if(!"R0" %in% names(dots)) stop("Sequence of R0 is needed for profile.")

  nll <- vapply(log(dots$R0 * Assessment@obj$env$data$rescale), function(xx) Assessment@obj$fn(xx), numeric(1))
  output <- new("prof", Model = Assessment@Model, Name = Assessment@Name, Par = "R0", MLE = Assessment@R0,
                grid = data.frame(R0 = dots$R0, nll = nll))
  return(output)
}


retrospective_SSS <- function(Assessment, nyr) {
  assign_Assessment_slots(Assessment)
  n_y <- info$data$n_y

  Year <- c(info$Year, max(info$Year) + 1)

  # Array dimension: Retroyr, Year, ts
  # ts includes: Calendar U, U_MSY, B, B/BMSY, B/B0, R, VB
  retro_ts <- array(NA, dim = c(nyr+1, n_y + 1, 7))
  TS_var <- c("U", "U_UMSY", "SSB", "SSB_SSBMSY", "SSB_SSB0", "R", "VB")
  dimnames(retro_ts) <- list(Peel = 0:nyr, Year = Year, Var = TS_var)

  SD <- summary(SD)
  retro_est <- array(NA, dim = c(nyr+1, dim(SD)))
  dimnames(retro_est) <- list(Peel = 0:nyr, Var = rownames(SD), Value = c("Estimate", "Std. Error"))

  lapply_fn <- function(i, info, obj) {
    n_y_ret <- n_y - i
    info$data$n_y <- n_y_ret
    info$data$C_hist <- info$data$C_hist[1:n_y_ret]

    dep <- info$data$I_hist[n_y]
    info$data$I_hist <- info$data$I_hist[1:n_y_ret]
    info$data$I_hist[n_y_ret] <- dep
    info$data$CAA_hist <- info$data$CAA_hist[1:n_y_ret, ]
    info$data$CAA_n <- info$data$CAA_n[1:n_y_ret]
    info$data$est_rec_dev <- info$data$est_rec_dev[1:n_y_ret]

    obj2 <- MakeADFun(data = info$data, parameters = info$params, map = obj$env$map,
                      DLL = "MSEtool", silent = TRUE)
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

      retro_ts[i+1, , ] <<- cbind(U, U_UMSY, SSB, SSB_SSBMSY, SSB_SSB0, R, VB)
      retro_est[i+1, , ] <<- summary(SD)

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

