
summary_SCA <- function(Assessment, SCA2 = FALSE) {
  assign_Assessment_slots(Assessment)

  if(conv) current_status <- c(F_FMSY[length(F_FMSY)], SSB_SSBMSY[length(SSB_SSBMSY)], SSB_SSB0[length(SSB_SSB0)])
  else current_status <- c(NA, NA, SSB_SSB0[length(SSB_SSB0)])
  current_status <- data.frame(Value = current_status)
  rownames(current_status) <- c("F/FMSY", "SSB/SSBMSY", "SSB/SSB0")

  Value <- c(h, info$data$M[1], info$data$max_age, info$LH$Linf, info$LH$K, info$LH$t0,
             info$LH$a * info$LH$Linf ^ info$LH$b, info$LH$A50, info$LH$A95)
  Description = c("Stock-recruit steepness", "Natural mortality", "Maximum age (plus-group)", "Asymptotic length", "Growth coefficient",
                  "Age at length-zero", "Asymptotic weight", "Age of 50% maturity", "Age of 95% maturity")
  rownam <- c("h", "M", "maxage", "Linf", "K", "t0", "Winf", "A50", "A95")
  input_parameters <- data.frame(Value = Value, Description = Description, stringsAsFactors = FALSE)
  rownames(input_parameters) <- rownam
  if(!"transformed_h" %in% names(obj$env$map) || (SCA2 && !info$fix_h)) input_parameters <- input_parameters[-1, ]

  if(conv) Value <- c(VB0, SSB0, MSY, FMSY, VBMSY, SSBMSY, SSBMSY/SSB0)
  else Value <- rep(NA, 7)

  Description <- c("Unfished vulnerable biomass",
                   "Unfished spawning stock biomass (SSB)", "Maximum sustainable yield (MSY)", "Fishing mortality at MSY",
                   "Vulnerable biomass at MSY", "SSB at MSY", "Spawning depletion at MSY")
  derived <- data.frame(Value = Value, Description = Description, stringsAsFactors = FALSE)
  rownames(derived) <- c("VB0", "SSB0", "MSY", "FMSY", "VBMSY", "SSBMSY", "SSBMSY/SSB0")

  if(conv && SCA2 && !info$fix_h) {
    derived <- rbind(derived, c(h, "Stock-recruit steepness"))
    rownames(derived)[8] <- "h"
  }

  model_estimates <- sdreport_int(SD)
  if(!is.character(model_estimates)) {
    rownames(model_estimates)[rownames(model_estimates) == "logF"] <- paste0("logF_", names(FMort))
    rownames(model_estimates)[rownames(model_estimates) == "log_rec_dev"] <- paste0("log_rec_dev_", names(FMort)[as.logical(obj$env$data$est_rec_dev)])
  }

  output <- list(model = paste("Statistical Catch-at-Age", ifelse(SCA2, "(SCA2)", "(SCA)")),
                 current_status = current_status, input_parameters = input_parameters,
                 derived_quantities = derived, model_estimates = model_estimates,
                 log_likelihood = matrix(NLL, ncol = 1, dimnames = list(names(NLL), "Neg.LL")))
  return(output)
}

rmd_SCA <- function(Assessment, SCA2 = FALSE, ...) {
  ss <- rmd_summary(paste("Statistical Catch-at-Age", ifelse(SCA2, "(SCA2)", "(SCA)")))

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
  if(SCA2) {
    lead_par <- rmd_meanR(header = "## Assessment {.tabset}\n### Estimates and Model Fit\n")
  } else {
    lead_par <- c(rmd_R0(header = "## Assessment {.tabset}\n### Estimates and Model Fit\n"), rmd_h())
  }

  assess_fit <- c(lead_par,
                  rmd_sel(age, Assessment@Selectivity[nrow(Assessment@Selectivity), ], fig.cap = "Estimated selectivity at age."),
                  rmd_assess_fit("Index", "index"), rmd_assess_resid("Index"), rmd_assess_qq("Index", "index"),
                  rmd_assess_fit("Catch", "catch"), rmd_assess_resid("Catch"), rmd_assess_qq("Catch", "catch"),
                  rmd_fit_age_comps("bubble"), rmd_fit_age_comps("annual"),
                  rmd_residual("Dev", fig.cap = "Time series of recruitment deviations.", label = Assessment@Dev_type,
                               blue = any(as.numeric(names(Assessment@Dev)) < Assessment@info$Year[1])),
                  rmd_residual("Dev", "SE_Dev", fig.cap = "Time series of recruitment deviations with 95% confidence intervals.",
                               label = Assessment@Dev_type, conv_check = TRUE, blue = any(as.numeric(names(Assessment@Dev)) < Assessment@info$Year[1])))

  #### Time Series
  ts_output <- c(rmd_F(header = "### Time Series Output\n"), rmd_F_FMSY(), rmd_SSB(), rmd_SSB_SSBMSY(),
                 rmd_SSB_SSB0(), rmd_Kobe("SSB_SSBMSY", xlab = "expression(SSB/SSB[MSY])"), rmd_R(),
                 rmd_N(), rmd_N_at_age(), rmd_C_at_age(), rmd_C_mean_age())

  # Productivity
  Arec <- Assessment@TMB_report$Arec
  Brec <- Assessment@TMB_report$Brec
  SSB <- Assessment@SSB[1:(length(Assessment@SSB)-1)]

  SR <- ifelse(SCA2, Assessment@info$SR, Assessment@info$data$SR_type)
  if(SR == "BH") expectedR <- Arec * SSB / (1 + Brec * SSB) else {
    expectedR <- Arec * SSB * exp(-Brec * SSB)
  }
  estR <- Assessment@R[as.numeric(names(Assessment@R)) > Assessment@info$Year[1]]

  productivity <- c(rmd_SR(SSB, expectedR, estR, header = "### Productivity\n\n\n"),
                    rmd_SR(SSB, expectedR, estR, fig.cap = "Stock-recruit relationship (trajectory plot).", trajectory = TRUE),
                    rmd_yield_F("SCA"), rmd_yield_depletion("SCA"), rmd_sp())

  return(c(ss, LH_section, data_section, assess_fit, ts_output, productivity))
}


profile_likelihood_SCA <- function(Assessment, ...) {
  dots <- list(...)
  if(!"R0" %in% names(dots) && !"h" %in% names(dots)) stop("Sequence of neither R0 nor h was found. See help file.")
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
    params$R0x <- log(profile_grid[i, 1]  * Assessment@obj$env$data$rescale)
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
                      inner.control = Assessment@info$inner.control, DLL = "MSEtool", silent = TRUE)
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


retrospective_SCA <- function(Assessment, nyr, SCA2 = FALSE) {
  assign_Assessment_slots(Assessment)
  n_y <- info$data$n_y

  Year <- c(info$Year, max(info$Year) + 1)

  # Array dimension: Retroyr, Year, ts
  # ts includes: Calendar F, F_MSY, B, B/BMSY, B/B0, R, VB
  retro_ts <- array(NA, dim = c(nyr+1, n_y + 1, 7))
  TS_var <- c("F", "F_FMSY", "SSB", "SSB_SSBMSY", "SSB_SSB0", "R", "VB")
  dimnames(retro_ts) <- list(Peel = 0:nyr, Year = Year, Var = TS_var)

  SD_nondev <- summary(SD)[rownames(summary(SD)) != "log_rec_dev" & rownames(summary(SD)) != "log_early_rec_dev"
                           & rownames(summary(SD)) != "logF", ]
  retro_est <- array(NA, dim = c(nyr+1, dim(SD_nondev)))
  dimnames(retro_est) <- list(Peel = 0:nyr, Var = rownames(SD_nondev), Value = c("Estimate", "Std. Error"))

  lapply_fn <- function(i, info, obj) {
    n_y_ret <- n_y - i
    info$data$n_y <- n_y_ret

    if(info$data$yindF + 1 > n_y_ret) {
      info_old <- info
      info$data$yindF <- as.integer(0.5 * n_y_ret)
      info$params$logF[info$data$yindF + 1] <- info_old$params$logF[info_old$data$yindF + 1]
    }

    info$data$C_hist <- info$data$C_hist[1:n_y_ret]
    info$data$I_hist <- info$data$I_hist[1:n_y_ret]
    info$data$CAA_hist <- info$data$CAA_hist[1:n_y_ret, ]
    info$data$CAA_n <- info$data$CAA_n[1:n_y_ret]
    info$data$est_rec_dev <- info$data$est_rec_dev[1:n_y_ret]

    info$params$log_rec_dev <- rep(0, n_y_ret)
    info$params$logF <- info$params$logF[1:n_y_ret]

    map <- obj$env$map
    if(any(names(map) == "log_rec_dev")) {
      new_map <- as.numeric(map$log_rec_dev) - i
      map$log_rec_dev <- factor(new_map[new_map > 0])
    }
    if(any(names(map) == "logF")) map$logF <- map$logF[1:n_y_ret]

    obj2 <- MakeADFun(data = info$data, parameters = info$params, map = map, random = obj$env$random,
                      inner.control = info$inner.control, DLL = "MSEtool", silent = TRUE)
    mod <- optimize_TMB_model(obj2, info$control)
    opt2 <- mod[[1]]
    SD <- mod[[2]]

    if(!is.character(opt2) && !is.character(SD)) {
      report <- obj2$report(obj2$env$last.par.best)
      if(SCA2) {
        ref_pt <- SCA_refpt_calc(E = report$E[1:(length(report$E) - 1)], R = report$R[2:length(report$R)], M = info$data$M,
                                 weight = info$data$weight, mat = info$data$mat, vul = report$vul, SR = info$SR,
                                 fix_h = info$fix_h, h = h)
      } else {
        ref_pt <- SCA_MSY_calc(Arec = report$Arec, Brec = report$Brec, M = info$data$M, weight = info$data$weight, mat = info$data$mat,
                               vul = report$vul, SR = info$data$SR_type)
      }

      report <- c(report, ref_pt)

      FMort <- c(report$F, rep(NA, i + 1))
      F_FMSY <- FMort/report$FMSY
      SSB <- c(report$E, rep(NA, i))
      SSB_SSBMSY <- SSB/report$EMSY
      SSB_SSB0 <- SSB/report$E0
      R <- c(report$R, rep(NA, i))
      VB <- c(report$E, rep(NA, i))
      #log_rec_dev <- c(report$log_rec_dev, rep(NA, i + 1))

      retro_ts[i+1, , ] <<- cbind(FMort, F_FMSY, SSB, SSB_SSBMSY, SSB_SSB0, R, VB)
      retro_est[i+1, , ] <<- summary(SD)[rownames(summary(SD)) != "log_rec_dev" & rownames(summary(SD)) != "log_early_rec_dev" &
                                          rownames(summary(SD)) != "logF", ]

      return(SD$pdHess)
    }
    return(FALSE)
  }

  conv <- vapply(0:nyr, lapply_fn, logical(1), info = info, obj = obj)
  if(any(!conv)) warning("Peels that did not converge: ", paste0(which(!conv) - 1, collapse = " "))

  retro <- new("retro", Model = Assessment@Model, Name = Assessment@Name, TS_var = TS_var, TS = retro_ts,
               Est_var = dimnames(retro_est)[[2]], Est = retro_est)
  attr(retro, "TS_lab") <- c("Fishing mortality", expression(F/F[MSY]), "Spawning biomass", expression(SSB/SSB[MSY]), "Spawning depletion",
                             "Recruitment", "Vulnerable biomass")

  return(retro)
}


summary_SCA2 <- function(Assessment) summary_SCA(Assessment, TRUE)

rmd_SCA2 <- function(Assessment) rmd_SCA(Assessment, TRUE)


profile_likelihood_SCA2 <- function(Assessment, ...) {
  dots <- list(...)
  if(!"meanR" %in% names(dots)) stop("Sequence of meanR was not found. See help file.")
  meanR <- dots$meanR

  params <- Assessment@info$params
  map <- Assessment@obj$env$map
  map$meanRx <- factor(NA)

  profile_fn <- function(i, Assessment, params, map) {

    params$meanRx <- log(meanR[i] * Assessment@obj$env$data$rescale)
    if(length(Assessment@opt$par) == 1) {
      nll <- Assessment@obj$fn(params$meanRx)
    } else {

      obj2 <- MakeADFun(data = Assessment@info$data, parameters = params, map = map, random = Assessment@obj$env$random,
                        inner.control = Assessment@info$inner.control, DLL = "MSEtool", silent = TRUE)
      opt2 <- optimize_TMB_model(obj2, Assessment@info$control)[[1]]
      if(!is.character(opt2)) nll <- opt2$objective else nll <- NA

    }
    return(nll)
  }
  nll <- vapply(1:length(meanR), profile_fn, numeric(1), Assessment = Assessment, params = params, map = map) - Assessment@opt$objective
  profile_grid <- data.frame(meanR = meanR, nll = nll)

  pars <- c("meanR")
  MLE <- Assessment@SD$value["meanR"]

  output <- new("prof", Model = Assessment@Model, Name = Assessment@Name, Par = pars, MLE = MLE, grid = profile_grid)
  return(output)
}


retrospective_SCA2 <- function(Assessment, nyr) retrospective_SCA(Assessment, nyr, TRUE)


plot_yield_SCA <- function(data, report, fmsy, msy, xaxis = c("F", "Biomass", "Depletion")) {
  xaxis <- match.arg(xaxis)
  F.vector = seq(0, 2.5 * fmsy, length.out = 1e2)

  M <- data$M
  mat <- data$mat
  weight <- data$weight
  maxage <- data$max_age
  SR <- data$SR_type

  vul <- report$vul

  BMSY <- report$EMSY
  B0 <- report$E0

  Arec <- report$Arec
  Brec <- report$Brec

  EPR <- Req <- NA
  solveMSY <- function(logF) {
    Fmort <- exp(logF)
    surv <- exp(-vul * Fmort - M)
    NPR <- c(1, cumprod(surv[1:(maxage-1)]))
    NPR[maxage] <- NPR[maxage]/(1 - surv[maxage])
    EPR <<- sum(NPR * mat * weight)
    if(SR == "BH") Req <<- (Arec * EPR - 1)/(Brec * EPR)
    if(SR == "Ricker") Req <<- log(Arec * EPR)/(Brec * EPR)
    CPR <- vul * Fmort/(vul * Fmort + M) * NPR * (1 - exp(-vul * Fmort - M))
    Yield <- Req * sum(CPR * weight)
    return(-1 * Yield)
  }

  Biomass <- Yield <- R <- rep(NA, length(F.vector))
  for(i in 1:length(F.vector)) {
    Yield[i] <- -1 * solveMSY(log(F.vector[i]))
    R[i] <- Req
    Biomass[i] <- EPR * Req
  }

  ind <- R >= 0

  if(xaxis == "F") {
    plot(F.vector[ind], Yield[ind], typ = 'l', xlab = "Fishing Mortality",
         ylab = "Equilibrium yield")
    segments(x0 = fmsy, y0 = 0, y1 = msy, lty = 2)
    segments(x0 = 0, y0 = msy, x1 = fmsy, lty = 2)
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
    segments(x0 = 0, y0 = msy, x1 = BMSY/B0, lty = 2)
    abline(h = 0, col = 'grey')
  }
  invisible(data.frame(F = F.vector[ind], Yield = Yield[ind], B = Biomass[ind], B_B0 = Biomass[ind]/report$B0))
}


