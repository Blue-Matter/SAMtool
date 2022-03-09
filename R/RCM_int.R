
RCM_int <- function(OM, RCMdata, condition = c("catch", "catch2", "effort"), selectivity = "logistic", s_selectivity = "B", LWT = list(),
                    comp_like = c("multinomial", "lognormal"), prior = list(),
                    max_F = 3, cores = 1L, integrate = FALSE, mean_fit = FALSE, drop_nonconv = FALSE,
                    drop_highF = FALSE,
                    control = list(iter.max = 2e+05, eval.max = 4e+05), ...) {
  
  dots <- list(...) # can be vul_par, ivul_par, log_rec_dev, log_early_rec_dev, map_vul_par, map_ivul_par, map_log_rec_dev, map_log_early_rec_dev, rescale, plusgroup, resample, OMeff, fix_dome
  if(!is.null(dots$maxF)) max_F <- dots$maxF
  
  comp_like <- match.arg(comp_like)
  condition <- match.arg(condition)
  
  dat_update <- check_RCMdata(RCMdata, OM, condition)
  OM <- dat_update$OM
  RCMdata <- dat_update$RCMdata
  StockPars <- dat_update$StockPars
  FleetPars <- dat_update$FleetPars
  ObsPars <- dat_update$ObsPars
  
  nsim <- OM@nsim
  proyears <- OM@proyears
  maxage <- OM@maxage
  nyears <- RCMdata@Misc$nyears
  nfleet <- RCMdata@Misc$nfleet
  nsurvey <- RCMdata@Misc$nsurvey
  
  OM@maxF <- max_F
  message("OM@maxF updated to ", max_F, ".")
  
  # No comp data
  if(!any(RCMdata@CAA > 0, na.rm = TRUE) && !any(RCMdata@CAL > 0, na.rm = TRUE)) {
    fix_sel <- TRUE
    message("No fishery length or age compositions were provided. Selectivity is fixed to values from OM.")
  } else {
    fix_sel <- FALSE
  }
  
  # Selectivity
  if(length(selectivity) == 1) selectivity <- rep(selectivity, RCMdata@Misc$nsel_block)
  if(length(selectivity) < RCMdata@Misc$nsel_block) stop("selectivity vector should be of length ", RCMdata@Misc$nsel_block, ").", call. = FALSE)
  sel <- int_sel(selectivity)
  
  # Survey selectivity
  if(nsurvey > 0) {
    if(is.null(s_selectivity)) s_selectivity <- rep("B", nsurvey)
    if(length(s_selectivity) == 1) s_selectivity <- rep(s_selectivity, nsurvey)
    s_sel <- int_s_sel(s_selectivity, nfleet)
  } else {
    s_sel <- int_s_sel("B")
  }
  
  # Likelihood weights
  RCMdata@Misc$LWT <- make_LWT(LWT, nfleet, nsurvey)
  
  # SR
  message(ifelse(OM@SRrel == 1, "Beverton-Holt", "Ricker"), " stock-recruitment relationship used.")
  
  # Generate priors
  prior <- make_prior(prior, nsurvey, OM@SRrel, dots)
  
  # Test for identical sims
  par_identical_sims <- par_identical_sims_fn(StockPars, FleetPars, ObsPars, RCMdata, dots)
  
  # Fit model
  if(!is.null(dots$resample) && dots$resample) { # Re-sample covariance matrix
    
    if(!requireNamespace("mvtnorm", quietly = TRUE)) stop("Please install the mvtnorm package.", call. = FALSE)
    
    message("\nResample = TRUE. Running mean fit model first...")
    mean_fit_output <- RCM_est(RCMdata = RCMdata, selectivity = sel, s_selectivity = s_sel,
                               LWT = RCMdata@Misc$LWT, comp_like = comp_like, prior = prior, 
                               max_F = max_F, integrate = integrate, StockPars = StockPars, ObsPars = ObsPars,
                               FleetPars = FleetPars, mean_fit = TRUE, dots = dots)
    
    if(length(mean_fit_output) > 0 && !mean_fit_output$report$conv) {
      warning("Mean fit model did not appear to converge. Will not be able to sample the covariance matrix.")
      message("Model did not converge. Returning the mean-fit model for evaluation.")
      
      drop_nonconv <- TRUE
      
      samps <- mean_fit_output$obj$env$last.par.best %>% matrix(nrow = nsim, ncol = length(mean_fit_output$obj$par), byrow = TRUE)
    } else {
      message("Model converged. Sampling covariance matrix for nsim = ", nsim, " replicates...")
      if(!all(par_identical_sims)) {
        message("Note: not all ", nsim, " replicates are identical.")
        message("These parameters should be identical among simulations: ",
                which(!par_identical_sims) %>% names() %>% paste0(collapse = " "))
      }
      samps <- mvtnorm::rmvnorm(nsim, mean_fit_output$opt$par, mean_fit_output$SD$cov.fixed)
    }
    
    mod <- list(mean_fit_output)
    res <- lapply(1:nsim, RCM_report_samps, samps = samps, obj = mean_fit_output$obj, conv = mean_fit_output$report$conv)
    conv <- rep(mean_fit_output$report$conv, nsim)
    
  } else if(all(par_identical_sims)) { # All identical sims detected
      
    message("\nAll ", nsim, " replicates are identical. Fitting once and replicating single fit...")
    
    mean_fit_output <- RCM_est(RCMdata = RCMdata, selectivity = sel, s_selectivity = s_sel,
                               LWT = RCMdata@Misc$LWT, comp_like = comp_like, prior = prior, 
                               max_F = max_F, integrate = integrate, StockPars = StockPars, ObsPars = ObsPars,
                               FleetPars = FleetPars, mean_fit = TRUE, control = control, dots = dots)
    
    mod <- list(mean_fit_output)
    res <- lapply(1:nsim, function(x) mean_fit_output$report)
    conv <- rep(mean_fit_output$report$conv, nsim)
    
  } else {
    
    message("\nFitting model (", nsim, " simulations) ...")
    if(cores > 1 && !snowfall::sfIsRunning()) MSEtool::setup(as.integer(cores))
    if(snowfall::sfIsRunning()) {
      mod <- snowfall::sfClusterApplyLB(1:nsim, RCM_est, RCMdata = RCMdata, selectivity = sel, s_selectivity = s_sel,
                                        LWT = RCMdata@Misc$LWT, comp_like = comp_like, prior = prior, 
                                        max_F = max_F, integrate = integrate, StockPars = StockPars, ObsPars = ObsPars,
                                        FleetPars = FleetPars, control = control, dots = dots)
    } else {
      mod <- lapply(1:nsim, RCM_est, RCMdata = RCMdata, selectivity = sel, s_selectivity = s_sel,
                    LWT = RCMdata@Misc$LWT, comp_like = comp_like, prior = prior, 
                    max_F = max_F, integrate = integrate, StockPars = StockPars, ObsPars = ObsPars,
                    FleetPars = FleetPars, control = control, dots = dots)
    }
    
    if(mean_fit) { ### Fit to life history means if mean_fit = TRUE
      message("Generating additional model fit from mean values of parameters in the operating model...\n")
      mean_fit_output <- RCM_est(RCMdata = RCMdata, selectivity = sel, s_selectivity = s_sel,
                                 LWT = RCMdata@Misc$LWT, comp_like = comp_like, prior = prior, 
                                 max_F = max_F, integrate = integrate, StockPars = StockPars, ObsPars = ObsPars,
                                 FleetPars = FleetPars, mean_fit = TRUE, control = control, dots = dots)
      
      if(!mean_fit_output$report$conv) warning("Mean fit model did not appear to converge.")
    } else {
      mean_fit_output <- list()
    }
    
    res <- lapply(mod, getElement, "report")
    conv <- vapply(res, getElement, logical(1), name = "conv")
  }
  
  if(drop_highF) {
    highF <- vapply(res, function(x) max(getElement(x, "F"), na.rm = TRUE) >= max_F, logical(1))
    if(sum(highF)) message(sum(highF), " out of ", nsim , " model fits had F on the upper boundary (F = ", max_F, "; ", round(100*sum(highF)/nsim, 2), "% of simulations).\n")
    
    conv <- conv & !highF
  }
  if(drop_nonconv) {
    NaF <- vapply(res, function(x) any(is.na(x$F) | is.infinite(x$F)), logical(1))
    if(sum(NaF)) message(sum(NaF), " out of ", nsim , " iterations had F with NA's")
    
    conv <- conv & !NaF
  }
  if(sum(conv) < nsim) message("Non-converged iteration(s): ", paste(which(!conv), collapse = " "), "\n")
  if(!sum(conv)) {
    message("Non-converged for all iterations. Returning all for evaluation.")
    keep <- !logical(OM@nsim)
  } else if(sum(conv) < nsim && (drop_nonconv | drop_highF)) {
    message("Non-converged and/or highF iterations will be removed.\n")
    keep <- conv
  } else {
    keep <- !logical(OM@nsim)
  }
  
  ### Get parameters to update OM
  obj_data <- mod[[1]]$obj$env$data
  OM_par <- RCM_update_OM(res, obj_data, maxage, nyears, proyears)
  
  ### R0
  OM@cpars$R0 <- OM_par$R0
  message("Range of unfished age-0 recruitment (OM@cpars$R0): ", paste(round(range(OM@cpars$R0), 2), collapse = " - "))
  
  ### Depletion and init D - init D is only reported, OM setup for initD by adjusting rec devs
  message("Range of initial spawning depletion: ", paste(round(range(OM_par$initD), 2), collapse = " - "))
  
  OM@cpars$D <- OM_par$D
  message("Range of spawning depletion (OM@cpars$D): ", paste(round(range(OM@cpars$D), 2), collapse = " - "), "\n")
  
  ### Selectivity and F
  ### Find
  OM@isRel <- FALSE
  OM@cpars$V <- OM_par$V
  OM@cpars$Find <- OM_par$Find
  message("Historical F and selectivity trends set in OM@cpars$Find and OM@cpars$V, respectively.")
  message("Selectivity during projection period is set to that in most recent historical year.")
  
  OM@cpars$qs <- rep(1, nsim)
  Eff <- apply(OM@cpars$Find, 2, range)
  OM@EffLower <- Eff[1, ]
  OM@EffUpper <- Eff[2, ]
  if(length(OM@EffYears) != nyears) OM@EffYears <- 1:nyears
  if(length(OM@Esd) == 0 && is.null(OM@cpars$Esd)) OM@Esd <- c(0, 0)
  message("Historical effort trends set in OM@EffLower and OM@EffUpper.\n")
  
  ### Rec devs
  OM@cpars$Perr <- OM_par$procsd
  message("Recruitment standard deviation set in OM@cpars$Perr.")
  
  Perr_y <- StockPars$Perr_y
  Perr_y[, 1:maxage] <- OM_par$early_Perr
  Perr_y[, maxage + 1:nyears] <- OM_par$Perr
  message("Historical recruitment set in OM@cpars$Perr_y.")
  
  if(any(OM_par$AC != 0)) {
    OM@cpars$AC <- OM_par$AC
    OM@AC <- range(OM_par$AC)
    message("Range of recruitment autocorrelation OM@AC: ", paste(round(range(OM@AC), 2), collapse = " - "))
    
    OM@cpars$Perr_y <- RCM_sample_future_dev(obj_data$est_rec_dev, OM_par$procsd, OM_par$AC, 
                                             OM_par$log_rec_dev, Perr_y, maxage, nyears, proyears)
    message("Future recruitment deviations sampled with autocorrelation (in OM@cpars$Perr_y).\n")
  }
  
  ### Assign OM variables that were used in the RCM to the output
  OM@cpars$Len_age <- StockPars$Len_age
  OM@cpars$Linf <- StockPars$Linf
  OM@cpars$K <- StockPars$K
  OM@cpars$t0 <- StockPars$t0
  OM@cpars$LenCV <- StockPars$LenCV
  OM@cpars$LatASD <- StockPars$LatASD
  OM@cpars$Wt_age <- StockPars$Wt_age
  
  if(any(apply(StockPars$Mat_age, 1, function(x) all(x >= 0.5)))) { # Any simulations where all mat_age > 0.5?
    OM@cpars$L50 <- StockPars$L50
    OM@cpars$L95 <- StockPars$L95
  } else {
    OM@cpars$Mat_age <- StockPars$Mat_age
  }
  if(prior$use_prior[2]) {
    OM@cpars$hs <- OM_par$h
  } else {
    OM@cpars$hs <- StockPars$hs
  }
  if(prior$use_prior[3]) {
    OM@cpars$M_ageArray <- array(OM_par$Mest, c(nsim, maxage+1, nyears + proyears))
  } else {
    OM@cpars$M_ageArray <- StockPars$M_ageArray
  }
  
  if(any(RCMdata@CAL > 0, na.rm = TRUE) || (any(RCMdata@MS > 0, na.rm = TRUE) & RCMdata@MS_type == "length") ||
     any(RCMdata@IAL > 0, na.rm = TRUE)) {
    OM@cpars$CAL_bins <- RCMdata@Misc$lbin
    OM@cpars$CAL_binsmid <- RCMdata@Misc$lbinmid
    message("RCM length bins will be added to OM.")
  }
  
  if(!is.null(dots$plusgroup) && !dots$plusgroup) OM@cpars$plusgroup <- 0L
  if(!any(RCMdata@I_sd > 0, na.rm = TRUE)) OM@cpars$Iobs <- ObsPars$Iobs
  message("Growth, maturity, natural mortality, and steepness values from RCM are set in OM@cpars.\n")
  
  ### Output S4 object
  RCMdata@Misc$prior <- prior
  output <- new("RCModel", OM = MSEtool::SubCpars(OM, keep), 
                SSB = OM_par$SSB[keep, , drop = FALSE], 
                NAA = OM_par$NAA[keep, , , drop = FALSE], 
                CAA = OM_par$CAA[keep, , , , drop = FALSE], 
                CAL = OM_par$CAL[keep, , , , drop = FALSE], 
                mean_fit = mean_fit_output, 
                conv = conv[keep], 
                data = RCMdata, 
                Misc = res[keep])
  
  output@config <- list(drop_sim = which(!keep))
  
  # Data in cpars
  if(sum(RCMdata@Chist > 0, na.rm = TRUE) || nsurvey > 0) {
    real_Data <- new("Data")
    real_Data@Year <- (output@OM@CurrentYr - output@OM@nyears + 1):output@OM@CurrentYr
    real_Data@LHYear <- max(real_Data@Year)
    real_Data@MaxAge <- maxage
    if(sum(RCMdata@Chist > 0, na.rm = TRUE) && all(!is.na(RCMdata@Chist))) {
      real_Data@Cat <- matrix(rowSums(RCMdata@Chist, na.rm = TRUE), 1, nyears)
      real_Data@CV_Cat <- sqrt(exp((rowSums(RCMdata@C_sd * RCMdata@Chist)/rowSums(RCMdata@Chist))^2 - 1)) %>% 
        matrix(1, nyears)
      message("Historical catch data added to OM@cpars$Data@Cat.")
      if(nfleet > 1) message("Annual catch CV is a weighted average by fleet-specific catch.")
    }
    if(nfleet == 1) {
      if(sum(RCMdata@CAA, na.rm = TRUE)) {
        real_Data@CAA <- aperm(RCMdata@CAA, c(3, 1, 2))
        message("Age comps added to OM@cpars$Data@CAA (nfleet = 1).")
      }
      if(sum(RCMdata@CAL, na.rm = TRUE)) {
        real_Data@CAL <- aperm(RCMdata@CAL, c(3, 1, 2))
        real_Data@CAL_mids <- RCMdata@Misc$lbinmid
        real_Data@CAL_bins <- RCMdata@Misc$lbin
        message("Length comps added to OM@cpars$Data@CAL (nfleet = 1).")
      }
    }
    if(.hasSlot(real_Data, "AddInd") && nsurvey > 0) {
      real_Data@AddInd <- array(RCMdata@Index, c(nyears, nsurvey, output@OM@nsim)) %>%
        aperm(perm = c(3, 2, 1))
      real_Data@CV_AddInd <- array(sqrt(exp(RCMdata@I_sd^2) - 1), c(nyears, nsurvey, output@OM@nsim)) %>%
        aperm(perm = c(3, 2, 1))
      
      # Cannot accommodate indices mirrored to fleet when nfleet > 1 and fleet has time-varying sel
      real_Data@AddIndType <- vapply(s_sel, process_AddIndType, numeric(1), nfleet = nfleet)
      real_Data@AddIndV <- lapply(1:nsurvey, process_AddIndV, Misc = output@Misc, s_sel = s_sel,
                                  n_age = maxage + 1, nfleet = nfleet, nyears = nyears) %>%
        simplify2array() %>% aperm(c(1, 3, 2))
      real_Data@AddIunits <- RCMdata@I_units
      
      output@OM@cpars$AddIbeta <- matrix(1, output@OM@nsim, nsurvey)
      message("Historical indices added to OM@cpars$Data@AddInd.")
    }
    output@OM@cpars$Data <- real_Data
  }
  
  # Check whether observed matches predicted
  if(RCMdata@Misc$condition == "catch2") {
    catch_check_fn <- function(x, report, data) {
      if(report[[x]]$conv) {
        catch_diff <- report[[x]]$Cpred/RCMdata@Chist - 1
        flag <- max(abs(catch_diff), na.rm = TRUE) > 0.01 | any(is.na(report[[x]]$Cpred))
      } else {
        flag <- FALSE
      }
      return(flag)
    }
    
    do_catch_check <- vapply(1:length(res), catch_check_fn, logical(1), report = res, data = data)
    if(any(do_catch_check)) {
      flag_ind <- paste(which(do_catch_check), collapse = " ")
      message("Note: there is predicted catch that deviates from observed catch by more than 1% in simulations:")
      message(flag_ind)
    }
  }
  message("Complete.")
  
  return(output)
}


RCM_report_samps <- function(x, samps, obj, conv) {
  report <- obj$report(samps[x, ]) %>% RCM_posthoc_adjust(obj, par = samps[x, ])
  report$conv <- conv
  return(report)
}


RCM_update_OM <- function(res, obj_data, maxage, nyears, proyears) {
  out <- list()
  
  out$R0 <- vapply(res, getElement, numeric(1), "R0")
  out$initD <- vapply(res, function(x) x$E[1]/x$E0_SR, numeric(1))
  out$D <- vapply(res, function(x) x$E[length(x$E)-1]/x$E0_SR, numeric(1))
  
  make_F <- function(x) { # Extra step to avoid apical F = 0
    apicalF <- x$F
    apicalF[apicalF < 1e-4] <- 1e-4
    F_at_age <- lapply(1:ncol(apicalF), function(xx) apicalF[, xx] * x$vul[1:nyears, , xx]) %>% 
      simplify2array() %>% apply(1:2, sum)
    return(F_at_age)
  }
  F_matrix <- lapply(res, make_F)
  apical_F <- lapply(F_matrix, function(x) apply(x, 1, max))
  
  expand_V_matrix <- function(x) {
    y <- matrix(x[nyears, ], proyears, maxage + 1, byrow = TRUE)
    rbind(x, y)
  }
  
  out$V <- Map("/", e1 = F_matrix, e2 = apical_F) %>% lapply(expand_V_matrix) %>% simplify2array() %>% aperm(3:1)
  out$Find <- do.call(rbind, apical_F)
  
  out$procsd <- vapply(res, getElement, numeric(1), "tau")
  
  make_Perr <- function(x, obj_data) {
    bias_corr <- ifelse(obj_data$est_rec_dev, exp(-0.5 * x$tau^2), 1)
    res <- exp(x$log_rec_dev) * bias_corr
    res[1] <- res[1] * x$R_eq/x$R0
    return(res)
  }
  out$Perr <- sapply(res, make_Perr, obj_data = obj_data) %>% t()
  
  make_early_Perr <- function(x, obj_data) {
    res <- x$R_eq * x$NPR_equilibrium / x$R0 / x$NPR_unfished[1, ]
    bias_corr <- ifelse(obj_data$est_early_rec_dev, exp(-0.5 * x$tau^2), 1)
    early_dev <- exp(x$log_early_rec_dev) * bias_corr
    out <- res[-1] * early_dev
    return(rev(out))
  }
  out$early_Perr <- sapply(res, make_early_Perr, obj_data = obj_data) %>% t()
  
  out$log_rec_dev <- sapply(res, getElement, "log_rec_dev") %>% t()
  
  if(!all(out$log_rec_dev == 0)) {
    out$AC <- apply(out$log_rec_dev, 1, function(x) {
      out <- acf(x, lag.max = 1, plot = FALSE)$acf[2]
      ifelse(is.na(out), 0, out)
    })
  } else {
    out$AC <- rep(0, length(res))
  }
  
  out$h <- vapply(res, getElement, numeric(1), "h")
  out$Mest <- vapply(res, getElement, numeric(1), "Mest")
  
  out$SSB <- sapply(res, getElement, "E") %>% t()
  out$NAA <- sapply(res, getElement, "N", simplify = "array") %>% aperm(c(3, 1, 2))
  out$CAA <- sapply(res, getElement, "CAApred", simplify = "array") %>% aperm(c(4, 1:3))
  out$CAL <- sapply(res, getElement, "CALpred", simplify = "array") %>% aperm(c(4, 1:3))
  
  return(out)
}


RCM_sample_future_dev <- function(est_rec_dev, procsd, AC, log_rec_dev, Perr_y, maxage, nyears, proyears) {
  if(any(!est_rec_dev)) { # Sample historical rec devs for OM for most recent years
    yr_fixed_rec_dev <- which(!est_rec_dev)
    if(any(yr_fixed_rec_dev > sum(est_rec_dev))) {
      yr_hist_sample <- yr_fixed_rec_dev[yr_fixed_rec_dev > sum(est_rec_dev)] %>% min()
      nyears <- length(est_rec_dev)
      samp_hist <- Map(dev_AC, AC = AC, stdev = procsd, 
                       chain_start = log_rec_dev[, yr_hist_sample - 1],
                       MoreArgs = list(n = length(yr_hist_sample:nyears), mu = 1))
      log_rec_dev[, yr_hist_sample:nyears] <- do.call(rbind, samp_hist)
      
      Perr_y[, (maxage + yr_hist_sample):(maxage + nyears)] <- exp(log_rec_dev[, yr_hist_sample:nyears])
      message("Historical recruitment deviations sampled with autocorrelation starting in year ", yr_hist_sample, " out of OM@nyears = ", nyears)
    }
  }
  samp_proj <- Map(dev_AC, AC = AC, stdev = procsd, 
                   chain_start = log_rec_dev[, nyears],
                   MoreArgs = list(n = proyears, mu = 1))
  
  Perr_y[, maxage + nyears + 1:proyears] <- exp(do.call(rbind, samp_proj))
  return(Perr_y)
}

process_AddIndType <- function(s_sel, nfleet) {
  if(s_sel == -4 || s_sel == -2 || s_sel == -1 || s_sel == 0) { # -4 = B,-2 - 0 = custom sel
    return(1)
  } else if(s_sel == -3) { # SSB
    return(2)
  } else { # Fleet
    return(ifelse(nfleet > 1, 1, 3))
  }
}

process_AddIndV <- function(sur, Misc, s_sel, n_age, nfleet, nyears) { # Return a matrix of nsim x nages
  if(s_sel[sur] < -2 || (s_sel[sur] == 1 & nfleet == 1)) { # -4 = B, -3 = SSB, single-fleet VB
    out <- matrix(1, length(Misc), n_age)
  } else { # custom sel or multi-fleet
    out <- do.call(rbind, lapply(Misc, function(x) x$ivul[nyears, , sur]))
  }
  return(out)
}

RCM_retro <- function(x, nyr = 5) {
  if(length(x@mean_fit) == 0) stop("Re-run RCM() with argument `mean_fit = TRUE`", .call = FALSE)
  
  data <- x@mean_fit$obj$env$data
  params <- x@mean_fit$obj$env$parameters
  n_y <- data$n_y
  map <- x@mean_fit$obj$env$map
  
  if(data$nfleet > 1) {
    retro_ts <- array(NA, dim = c(nyr+1, n_y + 1, data$nfleet + 4))
    TS_var <- c(paste("Fleet", 1:data$nfleet, "F"), "Apical F", "SSB", "SSB_SSB0", "R")
  } else {
    retro_ts <- array(NA, dim = c(nyr+1, n_y + 1, 4))
    TS_var <- c("Apical F", "SSB", "SSB_SSB0", "R")
  }
  dimnames(retro_ts) <- list(Peel = 0:nyr, Year = (x@OM@CurrentYr - n_y):x@OM@CurrentYr + 1, Var = TS_var)
  
  new_args <- lapply(n_y - 0:nyr, RCM_retro_subset, data = data, params = params, map = map)
  lapply_fn <- function(i, new_args, x) {
    obj2 <- MakeADFun(data = new_args[[i+1]]$data, parameters = new_args[[i+1]]$params, map = new_args[[i+1]]$map,
                      random = x@mean_fit$obj$env$random, DLL = "SAMtool", silent = TRUE)
    if(new_args[[i+1]]$data$condition == "catch2") {
      if(any(is.na(obj2$report(obj2$par)$F)) || any(is.infinite(obj2$report(obj2$par)$F))) {
        for(ii in 1:10) {
          obj2$par["R0x"] <- 0.5 + obj2$par["R0x"]
          if(all(!is.na(obj2$report(obj2$par)$F)) && all(!is.infinite(obj2$report(obj2$par)$F))) break
        }
      }
    }
    mod <- optimize_TMB_model(obj2, control = list(iter.max = 2e+05, eval.max = 4e+05), restart = 0)
    opt2 <- mod[[1]]
    SD <- mod[[2]]
    
    if(!is.character(opt2)) {
      report <- obj2$report(obj2$env$last.par.best) %>% RCM_posthoc_adjust(obj2, dynamic_SSB0 = FALSE)
      
      FMort <- rbind(report$F, matrix(NA, i + 1, ncol(report$F)))
      if(data$nfleet > 1) FMort <- cbind(FMort, apply(report$F_at_age, 1, max, na.rm = TRUE) %>% c(rep(NA, i+1)))
      
      SSB <- c(report$E, rep(NA, i))
      SSB_SSB0 <- SSB/report$E0_SR
      R <- c(report$R, rep(NA, i))
      
      retro_ts[i+1, , ] <<- cbind(FMort, SSB, SSB_SSB0, R)
      
      return(SD$pdHess)
    }
    return(FALSE)
  }
  
  conv <- vapply(0:nyr, lapply_fn, logical(1), new_args = new_args, x = x)
  if(any(!conv)) warning("Peels that did not converge: ", paste0(which(!conv) - 1, collapse = " "))
  
  retro <- new("retro", Model = "RCM", Name = x@OM@Name, TS_var = TS_var, TS = retro_ts)
  if(data$nfleet > 1) {  
    attr(retro, "TS_lab") <- c(paste("Fishing mortality of Fleet", 1:data$nfleet), "Apical F",
                               "Spawning biomass", "Spawning depletion", "Recruitment")
  } else {  
    attr(retro, "TS_lab") <- c(paste("Fishing mortality of Fleet", 1:data$nfleet),
                               "Spawning biomass", "Spawning depletion", "Recruitment")
  }
  return(retro)
}


RCM_retro_subset <- function(yr, data, params, map) {
  data_out <- structure(data, check.passed = NULL)
  params_out <- lapply(params, function(x) if(!is.null(attr(x, "map"))) attr(x, "shape") else x)
  
  if(yr < data$n_y) {
    ##### Data object
    mat <- c("C_hist", "E_hist", "I_hist", "sigma_I", "sigma_C", "CAA_n", "CAL_n", "IAA_n", "IAL_n", "msize")
    if(!is.null(map$log_M)) mat <- c(mat, "M_data")
    mat_ind <- match(mat, names(data_out))
    data_out[mat_ind] <- lapply(data_out[mat_ind], function(x) x[1:yr, , drop = FALSE])
    
    mat2 <- c("len_age", "wt", "mat", "sel_block")
    mat_ind2 <- match(mat2, names(data_out))
    data_out[mat_ind2] <- lapply(data_out[mat_ind2], function(x) x[1:(yr+1), , drop = FALSE])
    
    # Update nsel_block, n_y
    data_out$nsel_block <- data_out$sel_block %>% as.vector() %>% unique() %>% length()
    data_out$n_y <- yr
    
    # Array 1:yr first index
    arr <- c("CAA_hist", "CAL_hist", "IAA_hist", "IAL_hist")
    arr_ind <- match(arr, names(data_out))
    data_out[arr_ind] <- lapply(data_out[arr_ind], function(x) x[1:yr, , , drop = FALSE])
    
    # Vector 1:yr
    data_out$est_rec_dev <- data_out$est_rec_dev[1:yr]
    
    ##### Parameters
    # Update rec devs
    params_out$log_rec_dev <- params_out$log_rec_dev[1:yr]
    if(!is.null(map$log_rec_dev)) map$log_rec_dev <- map$log_rec_dev[1:yr] %>% factor()
    
    # Update F
    if(data_out$condition == "catch") {
      params_out$log_F_dev <- params_out$log_F_dev[1:yr, , drop = FALSE]
      
      if(any(data_out$yind_F + 1 > yr)) {
        data_out$yind_F <- as.integer(0.5 * yr)
        params_out$log_F_dev[data_out$yind_F + 1, ] <- params$log_F_dev[data$yind_F + 1, ]
      }
    }
    
    ## Update nsel block if needed
    if(data$nsel_block > data_out$nsel_block) {
      sel_block_ind <- data_out$sel_block %>% as.vector() %>% unique()
      params_out$vul_par <- params_out$vul_par[, sel_block_ind, drop = FALSE]
      if(!is.null(map$vul_par)) map$vul_par <- map$vul_par[, sel_block_ind, drop = FALSE] %>% factor()
    }
  }
  
  return(list(data = data_out, params = params_out, map = map))
}


profile_likelihood_RCM <- function(x, ...) {
  dots <- list(...)
  if(!length(x@mean_fit)) stop("No model found. Re-run RCM with mean_fit = TRUE.")
  
  data <- x@mean_fit$obj$env$data
  params <- x@mean_fit$obj$env$parameters
  n_y <- data$n_y
  map <- x@mean_fit$obj$env$map
  LWT <- try(x@data@Misc$LWT, silent = TRUE)
  if(is.character(LWT)) LWT <- x@data$LWT
  
  new_args <- RCM_retro_subset(n_y, data = data, params = params, map = map)
  
  if(!is.null(dots$D)) {
    if(length(dots) > 1) message("Only doing depletion profile...")
    if(!requireNamespace("abind", quietly = TRUE)) {
      stop("Please install the abind package to run this profile.")
    }
    if(!requireNamespace("reshape2", quietly = TRUE)) {
      stop("Please install the reshape2 package to run this profile.")
    }
    profile_grid <- expand.grid(D = dots$D)
    
    new_args$data$IAA_hist <- abind::abind(new_args$data$IAA_hist, 
                                           array(NA, c(n_y, new_args$data$n_age, 1)), 
                                           along = 3)
    new_args$data$IAL_hist <- abind::abind(new_args$data$IAL_hist, 
                                           array(NA, c(n_y, length(new_args$data$lbinmid), 1)), 
                                           along = 3)
    new_args$data$IAA_n <- cbind(new_args$data$IAA_n, rep(0, n_y))
    new_args$data$IAL_n <- cbind(new_args$data$IAL_n, rep(0, n_y))
    new_args$data$nsurvey <- new_args$data$nsurvey + 1
    
    new_args$data$ivul_type <- c(new_args$data$ivul_type, -3)
    new_args$data$abs_I <- c(new_args$data$abs_I, 0)
    new_args$data$I_units <- c(new_args$data$I_units, 1) 
    new_args$data$LWT_index <- rbind(new_args$data$LWT_index, rep(1, 3))
    
    new_args$data$use_prior <- c(new_args$data$use_prior, 0)
    new_args$data$prior_dist <- rbind(new_args$data$prior_dist, rep(NA_real_, 2))
    
    new_args$params$ivul_par <- cbind(new_args$params$ivul_par, rep(0, nrow(new_args$params$ivul_par)))
    if(!is.null(new_args$map$ivul_par)) {
      new_args$map$ivul_par <- factor(c(new_args$map$ivul_par, rep(NA, nrow(new_args$params$ivul_par))))
    }
    
    LWT$Index <- c(LWT$Index, 1)
    LWT$IAA <- c(LWT$IAA, 1)
    LWT$IAL <- c(LWT$IAL, 1)
    
    profile_fn <- function(i, new_args, x) {
      # Add new survey of SSB depletion
      SSB_index <- SSB_Isd <- rep(NA_real_, n_y)
      SSB_index[1] <- 1
      SSB_index[n_y] <- profile_grid$D[i]
      SSB_Isd[c(1, n_y)] <- 0.01
      
      new_args$data$I_hist <- cbind(new_args$data$I_hist, SSB_index)
      new_args$data$sigma_I <- cbind(new_args$data$sigma_I, SSB_Isd)
      
      obj2 <- MakeADFun(data = new_args$data, parameters = new_args$params, map = new_args$map,
                        random = x@mean_fit$obj$env$random, DLL = "SAMtool", silent = TRUE)
      mod <- optimize_TMB_model(obj2, control = list(iter.max = 2e+05, eval.max = 4e+05), restart = 0)
      report <- obj2$report(obj2$env$last.par.best)
      RCM_get_likelihoods(report, LWT = LWT, f_name = paste0("Fleet_", 1:new_args$data$nfleet),
                          s_name = paste0("Index_", 1:(new_args$data$nsurvey-1)) %>% c("SSB_depletion"))
    }
    
    do_profile <- lapply(1:nrow(profile_grid), profile_fn, new_args = new_args, x = x)
    profile_grid$nll <- vapply(do_profile, function(xx) xx[[1]][1, 1], numeric(1))
    
    prof <- sapply(do_profile, nll_depletion_profile) %>% t()
    prof_out <- apply(prof, 2, sum) %>% as.logical()
    
    profile_grid <- cbind(profile_grid, prof[, prof_out] %>% as.data.frame())
    
    output <- new("prof", Model = "RCM", Par = "D", 
                  MLE = x@mean_fit$report$E[n_y]/x@mean_fit$report$E0_SR, grid = profile_grid)
  } else {
    
    if(is.null(dots$R0) && is.null(dots$h)) stop("Sequence of neither D, R0, nor h was found. See help file.")
    if(!is.null(dots$R0)) {
      R0 <- dots$R0 
    } else {
      R0 <- x@mean_fit$report$R0
      profile_par <- "h"
    }
    if(!is.null(dots$h)) {
      h <- dots$h 
    } else {
      h <- x@mean_fit$report$h
      profile_par <- "R0"
    }
    
    profile_grid <- expand.grid(R0 = R0, h = h)
    joint_profile <- !exists("profile_par", inherits = FALSE)
    
    profile_fn <- function(i, new_args, x) {
      new_args$params$R0x <- log(profile_grid$R0[i]  * new_args$data$rescale)
      if(new_args$data$SR_type == "BH") {
        new_args$params$transformed_h <- logit((profile_grid$h[i] - 0.2)/0.8)
      } else {
        new_args$params$transformed_h <- log(profile_grid$h[i] - 0.2)
      }
      
      if(joint_profile) {
        new_args$map$R0x <- new_args$map$transformed_h <- factor(NA)
      } else {
        if(profile_par == "R0") new_args$map$R0x <- factor(NA) else new_args$map$transformed_h <- factor(NA)
      }
      
      obj2 <- MakeADFun(data = new_args$data, parameters = new_args$params, map = new_args$map,
                        random = x@mean_fit$obj$env$random, DLL = "SAMtool", silent = TRUE)
      
      mod <- optimize_TMB_model(obj2, control = list(iter.max = 2e+05, eval.max = 4e+05), restart = 0)
      report <- obj2$report(obj2$env$last.par.best)
      RCM_get_likelihoods(report, LWT = LWT, f_name = paste0("Fleet_", 1:new_args$data$nfleet),
                          s_name = paste0("Index_", 1:new_args$data$nsurvey))
      
      #opt2 <- optimize_TMB_model(obj2, control = list(iter.max = 2e+05, eval.max = 4e+05), restart = 0)[[1]]
      #
      #if(!is.character(opt2)) nll <- opt2$objective else nll <- NA
      #return(nll)
    }
    
    do_profile <- lapply(1:nrow(profile_grid), profile_fn, new_args = new_args, x = x)
    profile_grid$nll <- vapply(do_profile, function(xx) xx[[1]][1, 1], numeric(1))
    
    prof <- sapply(do_profile, nll_depletion_profile) %>% t()
    prof_out <- apply(prof, 2, sum) %>% as.logical()
    
    profile_grid <- cbind(profile_grid, prof[, prof_out] %>% as.data.frame())
    
    #nll <- vapply(1:nrow(profile_grid), profile_fn, numeric(1), new_args = new_args, x = x)
    #profile_grid$nll <- nll
    
    if(joint_profile) {
      pars <- c("R0", "h")
      MLE <- c(x@mean_fit$report$R0, x@mean_fit$report$h)
    } else {
      pars <- profile_par
      MLE <- getElement(x@mean_fit$report, pars)
    }
    
    output <- new("prof", Model = "RCM", Par = pars, MLE = MLE, grid = profile_grid)
  }
  return(output)
}

#' @importFrom dplyr select starts_with
nll_depletion_profile <- function(xx) {
  nll_fleet <- xx[[2]][-nrow(xx[[2]]), ] %>% select(starts_with("Fleet"))
  nll_fleet$Data <- rownames(nll_fleet)
  nll_fleet <- reshape2::melt(nll_fleet, id.var = "Data", value.var = "nll")
  
  nll_index <- xx[[4]][-nrow(xx[[4]]), ] %>% select(starts_with("Index"))
  nll_index$Data <- rownames(nll_index)
  nll_index <- reshape2::melt(nll_index, id.var = "Data", value.var = "nll")
  
  nll <- rbind(nll_fleet, nll_index)
  
  nll$Name <- paste(nll$variable, nll$Data)
  
  c(nll$value, xx[[1]][c(2, 5, 6), 1]) %>% 
    structure(names = paste(nll$variable, nll$Data) %>% c("Recruitment Deviations", "Penalty", "Prior"))
}
