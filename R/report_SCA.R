
summary_SCA <- function(Assessment) {
  assign_Assessment_slots(Assessment)
  
  catch_eq <- obj$env$data$catch_eq
  SR <- obj$env$data$SR_type
  
  if (SR == "none" || !conv) {
    current_status <- data.frame(Value = c(ifelse(catch_eq == "Baranov", FMort[length(FMort)], U[length(U)]), 
                                           SSB[length(SSB)], SSB_SSB0[length(SSB_SSB0)]))
    rownames(current_status) <- c(ifelse(catch_eq == "Baranov", "F", "U"), "SSB", "SSB/SSB0")
  } else {
    current_status <- data.frame(Value = c(ifelse(catch_eq == "Baranov", F_FMSY[length(F_FMSY)], U_UMSY[length(U_UMSY)]), 
                                           SSB_SSBMSY[length(SSB_SSBMSY)], SSB_SSB0[length(SSB_SSB0)]))
    rownames(current_status) <- c(ifelse(catch_eq == "Baranov", "F/FMSY", "U/UMSY"), "SSB/SSBMSY", "SSB/SSB0")
  }
  
  if (SR == "none") h <- NA_real_
  Value <- c(h, TMB_report$M[1], info$data$n_age - 1, info$LH$Linf, info$LH$K, info$LH$t0,
             info$LH$a * info$LH$Linf ^ info$LH$b, info$LH$A50, info$LH$A95)
  Description = c("Stock-recruit steepness", "Natural mortality", "Maximum age (plus-group)", "Asymptotic length", "Growth coefficient",
                  "Age at length-zero", "Asymptotic weight", "Age of 50% maturity", "Age of 95% maturity")
  input_parameters <- data.frame(Value = Value, Description = Description, stringsAsFactors = FALSE)
  rownames(input_parameters) <- c("h", "M", "maxage", "Linf", "K", "t0", "Winf", "A50", "A95")
  if ("transformed_h" %in% names(obj$env$map) || SR != "none") {
    input_parameters <- input_parameters[-1, ]
  }

  if (conv && SR != "none") {
    Value <- c(VB0, SSB0, MSY, ifelse(catch_eq == "Pope", UMSY, FMSY), VBMSY, SSBMSY, SSBMSY/SSB0)
  } else {
    Value <- rep(NA_real_, 7)
  }
  Description <- c("Unfished vulnerable biomass",
                   "Unfished spawning stock biomass (SSB)", "Maximum sustainable yield (MSY)", 
                   ifelse(catch_eq == "Pope", "Exploitation rate at MSY", "Fishing mortality at MSY"),
                   "Vulnerable biomass at MSY", "SSB at MSY", "Spawning depletion at MSY")
  rownam <- c("VB0", "SSB0", "MSY", ifelse(catch_eq == "Pope", "UMSY", "FMSY"), "VBMSY", "SSBMSY", "SSBMSY/SSB0")
  if (conv && SR != "none" && "transformed_h" %in% names(obj$env$map)) {
    Description <- c(Description, "Stock-recruit steepness")
    Value <- c(Value, h)
    rownam <- c(rownam, "h")
  }
  derived <- data.frame(Value = Value, Description = Description, stringsAsFactors = FALSE)
  rownames(derived) <- rownam

  model_estimates <- sdreport_int(SD)
  if (!is.character(model_estimates)) {
    rownames(model_estimates)[rownames(model_estimates) %in% "logit_M_walk"] <- 
      paste0("logit_M_walk_", names(SSB)[-1])
    rownames(model_estimates)[rownames(model_estimates) %in% "logit_M"] <- 
      paste0("logit_M_", names(SSB))
    rownames(model_estimates)[rownames(model_estimates) %in% "log_F_dev"] <- 
      paste0("log_F_dev_", names(SSB)[-length(SSB)])
    rownames(model_estimates)[rownames(model_estimates) %in% "log_rec_dev"] <- 
      paste0("log_rec_dev_", names(FMort)[as.logical(obj$env$data$est_rec_dev)])
  }

  output <- list(model = "Statistical Catch-at-Age (SCA)",
                 current_status = current_status, input_parameters = input_parameters,
                 derived_quantities = derived, model_estimates = model_estimates,
                 log_likelihood = matrix(NLL, ncol = 1, dimnames = list(names(NLL), "Neg.LL")))
  return(output)
}

rmd_SCA <- function(Assessment, ...) {
  ss <- rmd_summary("Statistical Catch-at-Age (SCA)")
  
  age_comps <- any(Assessment@obj$env$data$CAA_n > 0)
  length_comps <- any(Assessment@obj$env$data$CAL_n > 0)

  # Life History
  LH_section <- c(rmd_LAA(header = "## Life History\n", SD_LAA = ifelse(length_comps, "info$LH$SD_LAA", "")), 
                  rmd_WAA(), rmd_LW(),
                  rmd_mat(fig.cap = "Maturity at age. Length-based maturity parameters were converted to the corresponding ages."))
  
  # Data section
  if (age_comps) {
    data_age_comps <- c(rmd_data_age_comps("bubble"), rmd_data_age_comps("annual"))
  } else {
    data_age_comps <- ""
  }
  if (length_comps) {
    data_length_comps <- c(rmd_data_length_comps("bubble"), rmd_data_length_comps("annual"))
  } else {
    data_length_comps <- ""
  }
  
  data_section <- c(rmd_data_timeseries("Catch", header = "## Data\n"), 
                    rmd_data_timeseries("Index", is_matrix = is.matrix(Assessment@Obs_Index), nsets = ncol(Assessment@Obs_Index)),
                    data_age_comps, data_length_comps)

  # Assessment
  #### Pars and Fit
  if (Assessment@obj$env$data$SR_type == "none") {
    lead_par <- rmd_R0(header = "## Assessment {.tabset}\n### Estimates and Model Fit\n")
  } else {
    lead_par <- c(rmd_R0(header = "## Assessment {.tabset}\n### Estimates and Model Fit\n"), rmd_h())
  }
  if (age_comps) {
    fit_age_comps <- c(rmd_fit_age_comps("bubble"), rmd_fit_age_comps("annual"))
  } else {
    fit_age_comps <- ""
  }
  if (length_comps) {
    fit_length_comps <- c(rmd_fit_length_comps("bubble"), rmd_fit_length_comps("annual"))
  } else {
    fit_length_comps <- ""
  }

  assess_fit <- c(lead_par, rmd_M_prior(), rmd_M_rw(),
                  rmd_sel(fig.cap = "Estimated selectivity at age."),
                  rmd_assess_fit("Catch", "catch"), rmd_assess_resid("Catch"), rmd_assess_qq("Catch", "catch"),
                  rmd_assess_fit_series(nsets = ncol(Assessment@Index)),
                  fit_age_comps, fit_length_comps,
                  rmd_residual("Dev", fig.cap = "Time series of recruitment deviations.", label = Assessment@Dev_type,
                               blue = any(as.numeric(names(Assessment@Dev)) < Assessment@info$Year[1])),
                  rmd_residual("Dev", "SE_Dev", fig.cap = "Time series of recruitment deviations with 95% confidence intervals.",
                               label = Assessment@Dev_type, conv_check = TRUE, 
                               blue = any(as.numeric(names(Assessment@Dev)) < Assessment@info$Year[1])))

  #### Time Series
  if (Assessment@obj$env$data$catch_eq == "Baranov") {
    F_output <- rmd_F(header = "### Time Series Output\n")
    if (Assessment@obj$env$data$SR_type != "none") F_output <- c(F_output, rmd_F_FMSY())
  } else {
    F_output <- rmd_U(header = "### Time Series Output\n")
    if (Assessment@obj$env$data$SR_type != "none") F_output <- c(F_output, rmd_U_UMSY())
  }
  
  if (age_comps) {
    C_age <- c(rmd_C_at_age(), rmd_C_mean_age())
  } else {
    C_age <- ""
  }
  if (length_comps) {
    C_length <- c(rmd_C_at_length(), rmd_C_mean_length())
  } else {
    C_length <- ""
  }
  
  ts_output <- c(F_output, rmd_M_rw(), rmd_M_DD(), rmd_SSB(),
                 rmd_dynamic_SSB0("TMB_report$dynamic_SSB0"), 
                 ifelse(Assessment@obj$env$data$SR_type != "none", rmd_SSB_SSBMSY(), ""),
                 rmd_SSB_SSB0(), 
                 ifelse(Assessment@obj$env$data$SR_type != "none", 
                        rmd_Kobe("SSB_SSBMSY", xlab = "expression(SSB/SSB[MSY])"), ""), 
                 rmd_R(), rmd_N(), rmd_N_at_age(), C_age, C_length)

  # Productivity
  if (Assessment@obj$env$data$SR_type != "none") {
    SR_calc <- c("SSB_SR <- SSB",
                 "if (info$data$SR_type == \"BH\") {",
                 "  R_SR <- TMB_report$Arec * SSB_SR / (1 + TMB_report$Brec * SSB_SR)",
                 "} else {",
                 "  R_SR <- TMB_report$Arec * SSB_SR * exp(-TMB_report$Brec * SSB_SR)",
                 "}",
                 "Rest <- R[as.numeric(names(R)) >= info$Year[1]]")
    SR_header <- c(rmd_SR(header = "### Productivity\n\n\n", SR_calc = SR_calc),
                   rmd_SR(fig.cap = "Stock-recruit relationship (trajectory plot).", trajectory = TRUE))
    if (Assessment@obj$env$data$catch_eq == "Baranov") {
      yield_curve <- c(rmd_yield_F("SCA"), rmd_yield_depletion("SCA"))
    } else {
      yield_curve <- c(rmd_yield_U("SCA_Pope"), rmd_yield_depletion("SCA_Pope"))
    }
    
  } else {
    SR_header <- "### Productivity\n\n\n"
    yield_curve <- ""
  }
  productivity <- c(SR_header, yield_curve, rmd_sp(yield_fn = Assessment@obj$env$data$SR_type != "none"), 
                    rmd_SPR(), rmd_YPR())

  return(c(ss, LH_section, data_section, assess_fit, ts_output, productivity))
}


profile_likelihood_SCA <- function(Assessment, ...) {
  dots <- list(...)
  if (!"R0" %in% names(dots) && !"h" %in% names(dots)) stop("Sequence of neither R0 nor h was found. See help file.")
  if (!is.null(dots$R0)) R0 <- dots$R0 else {
    R0 <- Assessment@R0
    profile_par <- "h"
  }
  if (!is.null(dots$h)) h <- dots$h else {
    h <- ifelse(length(Assessment@h), Assessment@h, NA_real_)
    profile_par <- "R0"
  }

  map <- Assessment@obj$env$map
  params <- Assessment@info$params

  profile_grid <- expand.grid(R0 = R0, h = h)
  joint_profile <- !exists("profile_par")

  profile_fn <- function(i, Assessment, params, map) {
    params$R0x <- log(profile_grid[i, 1]  * Assessment@obj$env$data$rescale)
    if (Assessment@info$data$SR_type == "BH") {
      params$transformed_h <- logit((profile_grid[i, 2] - 0.2)/0.8)
    } else if (Assessment@info$data$SR_type == "Ricker") {
      params$transformed_h <- log(profile_grid[i, 2] - 0.2)
    }

    if (joint_profile) {
      map$R0x <- map$transformed_h <- factor(NA)
    } else {
      if (profile_par == "R0") map$R0x <- factor(NA) else map$transformed_h <- factor(NA)
    }
    obj2 <- MakeADFun(data = Assessment@info$data, parameters = params, map = map, random = Assessment@obj$env$random,
                      inner.control = Assessment@info$inner.control, DLL = "SAMtool", silent = TRUE)
    opt2 <- optimize_TMB_model(obj2, Assessment@info$control, do_sd = FALSE)[[1]]
    if (!is.character(opt2)) nll <- opt2$objective else nll <- NA
    return(nll)
  }
  nll <- vapply(1:nrow(profile_grid), profile_fn, numeric(1), Assessment = Assessment, params = params, map = map) - Assessment@opt$objective
  profile_grid$nll <- nll

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


retrospective_SCA <- function(Assessment, nyr) { # Incorporates SCA, SCA2, and SCA_RWM
  assign_Assessment_slots(Assessment)
  n_y <- info$data$n_y
  
  Year <- c(info$Year, max(info$Year) + 1)
  
  # Array dimension: Retroyr, Year, ts
  # ts includes: Calendar F, B, R, VB
  if (any(grepl("logit_M", names(Assessment@SD$value)))) {
    TS_var <- c(ifelse(info$data$catch_eq == "Baranov", "F", "U"), "SSB", "R", "VB")
  } else {
    
    if (info$data$SR_type == "none") {
      TS_var <- c(ifelse(info$data$catch_eq == "Baranov", "F", "U"), "SSB", "SSB_SSB0", "R", "VB")
    } else {
      TS_var <- c(ifelse(info$data$catch_eq == "Baranov", "F", "U"), 
                  ifelse(info$data$catch_eq == "Baranov", "F_FMSY", "U_UMSY"), 
                  "SSB", "SSB_SSBMSY", "SSB_SSB0", "R", "VB")
    }
  }
  retro_ts <- array(NA, dim = c(nyr+1, n_y + 1, length(TS_var))) %>% 
    structure(dimnames = list(Peel = 0:nyr, Year = Year, Var = TS_var))
  
  SD_nondev <- summary(SD)[!rownames(summary(SD)) %in% 
                             c("log_rec_dev", "log_early_rec_dev", "log_F_dev", "logit_M", "logit_M_walk"), ]
  retro_est <- array(NA, dim = c(nyr+1, dim(SD_nondev))) %>% 
    structure(dimnames = list(Peel = 0:nyr, Var = rownames(SD_nondev), Value = c("Estimate", "Std. Error")))
  
  lapply_fn <- function(i, info, obj) {
    n_y_ret <- n_y - i
    info$data$n_y <- n_y_ret
    
    if (info$data$yindF + 1 > n_y_ret) {
      info_old <- info
      info$data$yindF <- as.integer(0.5 * n_y_ret)
      info$params$log_F_dev[info$data$yindF + 1] <- info_old$params$log_F_dev[info_old$data$yindF + 1]
    }
    
    info$data$C_hist <- info$data$C_hist[1:n_y_ret]
    
    if (Assessment@Model == "SSS") dep <- info$data$I_hist[n_y]
    info$data$I_hist <- info$data$I_hist[1:n_y_ret, , drop = FALSE]
    if (Assessment@Model == "SSS") info$data$I_hist[n_y_ret, ] <- dep
    
    info$data$CAA_hist <- info$data$CAA_hist[1:n_y_ret, ]
    info$data$CAA_n <- info$data$CAA_n[1:n_y_ret]
    
    info$data$CAL_hist <- info$data$CAL_hist[1:n_y_ret, , drop = FALSE]
    info$data$CAL_n <- info$data$CAL_n[1:n_y_ret]
    
    info$data$est_rec_dev <- info$data$est_rec_dev[1:n_y_ret]
    
    info$params$log_rec_dev <- rep(0, n_y_ret)
    info$params$log_F_dev <- info$params$log_F_dev[1:n_y_ret]
    info$params$logit_M_walk <- rep(0, n_y_ret)
    
    map <- obj$env$map
    if (any(names(map) == "log_rec_dev")) {
      new_map <- as.numeric(map$log_rec_dev) - i
      if (all(is.na(new_map))) {
        map$log_rec_dev <- factor(rep(NA, n_y_ret))
      } else {
        map$log_rec_dev <- factor(new_map[new_map > 0])
      }
    }
    if (any(names(map) == "log_F_dev")) map$log_F_dev <- map$log_F_dev[1:n_y_ret]
    if (any(names(map) == "logit_M_walk")) map$logit_M_walk <- map$logit_M_walk[1:n_y_ret]
    
    obj2 <- MakeADFun(data = info$data, parameters = info$params, map = map, random = obj$env$random,
                      inner.control = info$inner.control, DLL = "SAMtool", silent = TRUE)
    mod <- optimize_TMB_model(obj2, info$control)
    opt2 <- mod[[1]]
    SD <- mod[[2]]
    
    if (!is.character(opt2)) {
      report <- obj2$report(obj2$env$last.par.best)
      
      if (info$data$SR_type != "none") {
        ref_pt <- ref_pt_SCA(obj = obj2, report = report)
      }
      
      if ("F" %in% TS_var) retro_ts[i+1, , TS_var == "F"] <<- c(report$F, rep(NA, i + 1))
      if ("F_FMSY" %in% TS_var) retro_ts[i+1, , TS_var == "F_FMSY"] <<- c(report$F/ref_pt$FMSY, rep(NA, i + 1))
      if ("U" %in% TS_var) retro_ts[i+1, , TS_var == "U"] <<- c(report$U, rep(NA, i + 1))
      if ("U_UMSY" %in% TS_var) retro_ts[i+1, , TS_var == "U_UMSY"] <<- c(report$U/ref_pt$UMSY, rep(NA, i + 1))
      
      if ("SSB_SSBMSY" %in% TS_var) retro_ts[i+1, , TS_var == "SSB_SSBMSY"] <<- c(report$E/ref_pt$EMSY, rep(NA, i))
      if ("SSB_SSB0" %in% TS_var) retro_ts[i+1, , TS_var == "SSB_SSB0"] <<- c(report$E/report$E0, rep(NA, i))
      
      retro_ts[i+1, , TS_var == "SSB"] <<- c(report$E, rep(NA, i))
      retro_ts[i+1, , TS_var == "R"] <<- c(report$R, rep(NA, i))
      retro_ts[i+1, , TS_var == "VB"] <<- c(report$VB, rep(NA, i))
      
      retro_est[i+1, , ] <<- summary(SD)[!rownames(summary(SD)) %in% 
                                           c("log_rec_dev", "log_early_rec_dev", "log_F_dev", "logit_M", "logit_M_walk"), ]
      
      return(SD$pdHess)
    }
    return(FALSE)
  }
  
  conv <- vapply(0:nyr, lapply_fn, logical(1), info = info, obj = obj)
  if (any(!conv)) warning("Peels that did not converge: ", paste0(which(!conv) - 1, collapse = " "))
  
  retro <- new("retro", Model = Assessment@Model, Name = Assessment@Name, TS_var = TS_var, TS = retro_ts,
               Est_var = dimnames(retro_est)[[2]], Est = retro_est)
  
  TS_master_var <- c("F", "F_FMSY", "U", "U_UMSY", "SSB", "SSB_SSBMSY", "SSB_SSB0", "R", "VB")
  TS_master_lab <- c("Fishing mortality", expression(F/F[MSY]), "Exploitation rate", expression(U/U[MSY]),
                     "Spawning biomass", expression(SSB/SSB[MSY]), "Spawning depletion",
                     "Recruitment", "Vulnerable biomass")
  attr(retro, "TS_lab") <- TS_master_lab[match(TS_var, TS_master_var)]
  
  return(retro)
}

plot_yield_SCA <- function(data, report, fmsy, msy, xaxis = c("F", "Biomass", "Depletion")) {
  xaxis <- match.arg(xaxis)
  F.vector = seq(0, 2.5 * fmsy, length.out = 1e2)
  
  yield <- lapply(F.vector, yield_fn_SCA, M = report$M[nrow(report$M), ], mat = data$mat, weight = data$weight, vul = report$vul,
                  SR = data$SR_type, Arec = report$Arec, Brec = report$Brec, catch_eq = "Baranov", opt = FALSE,
                  B0 = report$B0, tv_M = data$tv_M, M_bounds = data$M_bounds)
  
  Biomass <- vapply(yield, getElement, numeric(1), "E")
  Yield <- vapply(yield, getElement, numeric(1), "Yield")
  R <- vapply(yield, getElement, numeric(1), "R")
  ind <- R >= 0
  
  BMSY <- report$EMSY
  B0 <- report$E0

  if (xaxis == "F") {
    plot(F.vector[ind], Yield[ind], typ = 'l', xlab = "Fishing Mortality",
         ylab = "Equilibrium yield")
    segments(x0 = fmsy, y0 = 0, y1 = msy, lty = 2)
    segments(x0 = 0, y0 = msy, x1 = fmsy, lty = 2)
    abline(h = 0, col = 'grey')
  }

  if (xaxis == "Biomass") {
    plot(Biomass[ind], Yield[ind], typ = 'l', xlab = "Spawning Stock Biomass",
         ylab = "Equilibrium yield")
    segments(x0 = BMSY, y0 = 0, y1 = msy, lty = 2)
    segments(x0 = 0, y0 = msy, x1 = BMSY, lty = 2)
    abline(h = 0, col = 'grey')
  }

  if (xaxis == "Depletion") {
    plot(Biomass[ind]/B0, Yield[ind], typ = 'l',
         xlab = expression(SSB/SSB[0]), ylab = "Equilibrium yield")
    segments(x0 = BMSY/B0, y0 = 0, y1 = msy, lty = 2)
    segments(x0 = 0, y0 = msy, x1 = BMSY/B0, lty = 2)
    abline(h = 0, col = 'grey')
  }
  invisible(data.frame(F = F.vector[ind], Yield = Yield[ind], B = Biomass[ind], B_B0 = Biomass[ind]/B0))
}


plot_yield_SCA_Pope <- function(data, report, umsy, msy, xaxis = c("U", "Biomass", "Depletion")) {
  xaxis <- match.arg(xaxis)
  if (xaxis == "U") {
    u.vector = seq(0, max(1, 2.5 * umsy), length.out = 100) 
  } else {
    u.vector = seq(0, 1, length.out = 100)
  }
  
  yield <- lapply(u.vector, yield_fn_SCA, M = report$M, mat = data$mat, weight = data$weight, vul = report$vul, 
                  SR = data$SR_type, Arec = report$Arec, Brec = report$Brec, catch_eq = "Pope", opt = FALSE,
                  B0 = report$B0, tv_M = data$tv_M, M_bounds = data$M_bounds)
  
  Biomass <- vapply(yield, getElement, numeric(1), "E")
  Yield <- vapply(yield, getElement, numeric(1), "Yield")
  R <- vapply(yield, getElement, numeric(1), "R")
  ind <- R >= 0
  
  BMSY <- report$EMSY
  B0 <- report$E0
  
  if (xaxis == "U") {
    plot(u.vector[ind], Yield[ind], typ = 'l', xlab = "Exploitation rate (U)",
         ylab = "Equilibrium yield")
    segments(x0 = umsy, y0 = 0, y1 = msy, lty = 2)
    segments(x0 = 0, y0 = msy, x1 = umsy, lty = 2)
    abline(h = 0, col = 'grey')
  }
  
  if (xaxis == "Biomass") {
    plot(Biomass[ind], Yield[ind], typ = 'l', xlab = "Spawning Stock Biomass",
         ylab = "Equilibrium yield")
    segments(x0 = BMSY, y0 = 0, y1 = msy, lty = 2)
    segments(x0 = 0, y0 = msy, x1 = BMSY, lty = 2)
    abline(h = 0, col = 'grey')
  }
  
  if (xaxis == "Depletion") {
    plot(Biomass[ind]/B0, Yield[ind], typ = 'l',
         xlab = expression(SSB/SSB[0]), ylab = "Equilibrium yield")
    segments(x0 = BMSY/B0, y0 = 0, y1 = msy, lty = 2)
    segments(x0 = 0, y0 = msy, x1 = BMSY/B0, lty = 2)
    abline(h = 0, col = 'grey')
  }
  invisible(data.frame(U = u.vector[ind], Yield = Yield[ind], B = Biomass[ind], B_B0 = Biomass[ind]/B0))
}
